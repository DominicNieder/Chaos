// Slide 3: the solid energy shell cut by a moving plane x = c. The cut face is the
// exact cross-section C(c) = { ½p_y² + V(c, y) ≤ E } and carries ALL crossings of all
// trajectories with the plane (filled: p_x > 0, hollow: p_x ≤ 0). No trajectory
// lines by default ("T" shows them briefly, with the shell as a cutaway).
//
// Slide 4: Poincaré map at x = 0 (unchanged 2D canvas, now rendered sharp at any size).

import { drawSectionMapAxes } from "./axis-indicator.js";
import {
  DEBUG, THREE, animate, attachRotation, clamp, createLineMaterial, createPointMaterial,
  createPolyline, createShellMaterial, createTrajectoryLine, createView, crossingsAt, fitPointsFromGeometry,
  handover, insetsFor, loadGroup, loadManifest, onSlide, registerSlideKey, revealScale,
  shellGeometry, whenReveal,
} from "./hh-scene.js";
import { buildCapGeometry, shellBounds, shellXExtent } from "./hh-energy.js";

// ---- Settings (tune on the real projector) -------------------------------------------
const SETTINGS = {
  keepSide: "auto",            // "auto": keep the half behind the cut face; or "above" (x ≥ c) / "below" (x ≤ c)
  chaoticPointSize: 3.4,       // CSS px of the slide
  regularPointSize: 6.5,
  chaoticBrightness: 0.85,
  regularBrightness: 1.0,
  hollowFadedOpacity: 0,       // opacity of p_x ≤ 0 points after "D" (teaser for slide 4)
  pointOffset: 0.004,          // world units towards the viewer, keeps points above the cap
  defaultCamera: { azimuth: 1.2, elevation: 0.9, zoom: 0.95 },   // looks onto the cut face (former default 0.9 / 0.48)
  entryAnimationMs: 1200,
  shellColor: "#3d4b5e",
  capColor: "#5f7088",
  capOutlineColor: "#e8ebf0",
  planeFrameColor: "#717a89",
  lineWidth: 1.4,
  planeRange: null,
  planeStep: 0.001,
  planeStops: [0],
  planeStart: "edge"
};

// ===================================================================================
// Slide 3 – moving section through the energy shell
// ===================================================================================
/** Add the 8 corners of the plane frame's sweep so the frame stays in view. */
function withFrameCorners(points, xExtent, bounds) {
  const corners = [];
  for (const x of [-xExtent, xExtent]) {
    for (const y of bounds.y) for (const py of bounds.py) corners.push(x, y, py);
  }
  const result = new Float32Array(points.length + corners.length);
  result.set(points);
  result.set(corners, points.length);
  return result;
}

/** Legend from the manifest: each chaotic orbit by name, regular orbits grouped. */
function buildOrbitLegend(target, entries) {
  if (!target) return;
  target.replaceChildren();
  const swatch = color => {
    const element = document.createElement("i");
    element.className = "hh-legend-swatch";
    element.style.background = color;
    return element;
  };
  const addGroup = (colors, label) => {
    const span = document.createElement("span");
    span.style.display = "inline-flex";
    span.style.alignItems = "center";
    span.style.gap = "0.2em";
    span.style.marginRight = "0.9em";
    colors.forEach(color => span.append(swatch(color)));
    span.append(document.createTextNode(` ${label}`));
    target.append(span);
  };
  for (const entry of entries.filter(e => e.kind === "chaotic")) addGroup([entry.color], entry.label ?? "Chaotic");
  const regular = entries.filter(e => e.kind !== "chaotic");
  if (regular.length) addGroup(regular.map(e => e.color), "Regular");
}

async function initializeMovingSection() {
  const canvas = document.querySelector("#hh-moving-section-plot");
  const slider = document.querySelector("#hh-moving-section-position");
  const valueLabel = document.querySelector("#hh-moving-section-value");
  const status = document.querySelector("#hh-moving-section-status");
  if (!canvas || !slider || !valueLabel || !status) return;
  const container = canvas.closest(".hh-viz");
  const legend = container.querySelector(".hh-viz-legend");
  const hollowButton = container.querySelector(".hh-toggle-hollow");

  try {
    const manifest = await loadManifest();
    const E = manifest.energy;
    const entries = await loadGroup("poincare-map");
    const xExtent = shellXExtent(E);
    const step = SETTINGS.planeStep;
    const snap = value => (step ? Math.round(value / step) * step : value);
    let [xMin, xMax] = SETTINGS.planeRange
      ? [Math.max(SETTINGS.planeRange[0], -xExtent), Math.min(SETTINGS.planeRange[1], xExtent)]
      : [-xExtent, xExtent];
    if (step) {                       // limits on the step grid, so the grid contains 0
      xMin = Math.ceil(xMin / step - 1e-9) * step;
      xMax = Math.floor(xMax / step + 1e-9) * step;
    }
    const bounds = shellBounds(E);
    const shell = shellGeometry(E);
    buildOrbitLegend(legend?.querySelector(".hh-legend-orbits"), entries);

    const camera = { ...SETTINGS.defaultCamera };
    const view = createView({
      container,
      canvas,
      axisCanvas: container.querySelector(".hh-axis-overlay"),
      camera,
      fitPoints: fitPointsFromGeometry(shell),
      insets: insetsFor(container),
    });

    // Clipping plane: keeps x ≥ c ("above") or x ≤ c ("below").
    const clipPlane = new THREE.Plane(new THREE.Vector3(1, 0, 0), 0);
    let keepBelow = false;
    let lockedSide = null;           // fixed side during the entry animation
    let c = xMin;

    function updateClipPlane() {
      if (keepBelow) clipPlane.set(new THREE.Vector3(-1, 0, 0), c);
      else clipPlane.set(new THREE.Vector3(1, 0, 0), -c);
    }

    // Solid shell (clipped) -----------------------------------------------------------
    const shellMaterial = createShellMaterial({ color: SETTINGS.shellColor, clippingPlanes: [clipPlane] });
    const shellMesh = view.add(new THREE.Mesh(shell, shellMaterial));

    // Cut face (exact cross-section), rebuilt when c changes ----------------------------
    const CAP_SAMPLES = 160;
    const capGeometry = new THREE.BufferGeometry();
    const capPositions = new Float32Array((CAP_SAMPLES + 1) * 2 * 3);
    const capNormals = new Float32Array((CAP_SAMPLES + 1) * 2 * 3);
    for (let i = 0; i < (CAP_SAMPLES + 1) * 2; i += 1) capNormals[3 * i] = 1;
    capGeometry.setAttribute("position", new THREE.BufferAttribute(capPositions, 3));
    capGeometry.setAttribute("normal", new THREE.BufferAttribute(capNormals, 3));
    const capMaterial = new THREE.MeshStandardMaterial({
      color: SETTINGS.capColor, roughness: 0.85, metalness: 0, side: THREE.DoubleSide,
      polygonOffset: true, polygonOffsetFactor: 1, polygonOffsetUnits: 2,
    });
    const capMesh = view.add(new THREE.Mesh(capGeometry, capMaterial));
    capMesh.frustumCulled = false;
    const outlineMaterial = view.trackLineMaterial(createLineMaterial({
      color: SETTINGS.capOutlineColor, linewidth: 1.6, depthMin: 1,
    }));
    const capOutline = view.add(createPolyline(new Float32Array(6), outlineMaterial));
    capOutline.renderOrder = 6;
    const frameMaterial = view.trackLineMaterial(createLineMaterial({
      color: SETTINGS.planeFrameColor, linewidth: 1, depthMin: 1,
    }));
    const planeFrame = view.add(createPolyline(new Float32Array(6), frameMaterial));

    // Crossing points ------------------------------------------------------------------
    // Draw order: chaotic below regular, hollow below filled.
    const groups = entries.map(entry => {
      const capacity = Math.max(1, entry.file.index.maxBin);
      const chaotic = entry.kind === "chaotic";
      const make = (style, order) => {
        const geometry = new THREE.BufferGeometry();
        geometry.setAttribute("position", new THREE.BufferAttribute(new Float32Array(3 * capacity), 3));
        geometry.setDrawRange(0, 0);
        const material = view.trackPointMaterial(createPointMaterial({
          color: entry.color,
          size: chaotic ? SETTINGS.chaoticPointSize : SETTINGS.regularPointSize,
          brightness: chaotic ? SETTINGS.chaoticBrightness : SETTINGS.regularBrightness,
          style,
          offset: SETTINGS.pointOffset,
        }));
        const points = view.add(new THREE.Points(geometry, material));
        points.frustumCulled = false;
        points.renderOrder = order;
        return points;
      };
      return {
        entry,
        buffer: new Float32Array(3 * capacity),
        filled: make("filled", chaotic ? 3 : 5),
        hollow: make("hollow", chaotic ? 2 : 4),
        line: null,
      };
    });

    function updateCrossings() {
      const start = DEBUG ? performance.now() : 0;
      let total = 0;
      for (const group of groups) {
        const n = crossingsAt(c, group.entry.file, group.buffer);
        const filled = group.filled.geometry.attributes.position;
        const hollow = group.hollow.geometry.attributes.position;
        let nf = 0;
        let nh = 0;
        for (let i = 0; i < n; i += 1) {
          const y = group.buffer[3 * i];
          const py = group.buffer[3 * i + 1];
          if (group.buffer[3 * i + 2] > 0) { filled.array.set([c, y, py], 3 * nf); nf += 1; }
          else { hollow.array.set([c, y, py], 3 * nh); nh += 1; }
        }
        filled.needsUpdate = true;
        hollow.needsUpdate = true;
        filled.clearUpdateRanges();
        hollow.clearUpdateRanges();
        filled.addUpdateRange(0, 3 * nf);
        hollow.addUpdateRange(0, 3 * nh);
        group.filled.geometry.setDrawRange(0, nf);
        group.hollow.geometry.setDrawRange(0, nh);
        total += n;
      }
      if (DEBUG) console.debug(`[hh] crossings at x = ${c.toFixed(4)}: ${total} in ${(performance.now() - start).toFixed(2)} ms`);
    }

    function updateCap() {
      const cap = buildCapGeometry(c, E, CAP_SAMPLES);
      const hasCap = Boolean(cap.interval);
      capMesh.visible = hasCap && cap.indices.length > 0 && !hint;
      capOutline.visible = hasCap;
      if (hasCap) {
        capPositions.set(cap.positions);
        capGeometry.attributes.position.needsUpdate = true;
        if (cap.indices.length && !capGeometry.index) {
          capGeometry.setIndex(new THREE.BufferAttribute(cap.indices, 1));
        }
        const outline = cap.outline.length === 3
          ? new Float32Array([...cap.outline, ...cap.outline])
          : cap.outline;
        capOutline.setPositions(outline);
      }
      const [y0, y1] = bounds.y;
      const [p0, p1] = bounds.py;
      planeFrame.setPositions(new Float32Array([c, y0, p0, c, y1, p0, c, y1, p1, c, y0, p1, c, y0, p0]));
    }

    function setPlane(value, { jump = false } = {}) {
      let next = clamp(snap(value), xMin, xMax);
      // Never jump over a stop (e.g. x = 0) while the user moves the plane.
      if (!jump) {
        for (const stop of SETTINGS.planeStops ?? []) {
          if ((c - stop) * (next - stop) < 0) { next = stop; break; }
        }
      }
      c = Object.is(next, -0) ? 0 : next;
      slider.value = String(c);
      valueLabel.textContent = `x = ${Number(c.toPrecision(3))}`;
      updateClipPlane();
      updateCap();
      updateCrossings();
      view.invalidate();
    }

    view.onBeforeRender(basis => {
      const side = lockedSide ?? SETTINGS.keepSide;
      const below = side === "auto" ? basis.e3.x >= 0 : side === "below";
      if (below !== keepBelow) { keepBelow = below; updateClipPlane(); }
    });

    // ---- Optional context: trajectory lines inside a cutaway shell ("T") -------------
    let hint = false;
    function setHint(value) {
      hint = value;
      for (const group of groups) {
        if (hint && !group.line) {
          const material = view.trackLineMaterial(createLineMaterial({
            color: group.entry.color, linewidth: SETTINGS.lineWidth, depthMin: 0.4, clippingPlanes: [clipPlane],
          }));
          group.line = view.add(createTrajectoryLine(group.entry.file, material));
          group.line.renderOrder = 1;
        }
        if (group.line) group.line.visible = hint;
      }
      shellMaterial.side = hint ? THREE.BackSide : THREE.FrontSide;
      shellMaterial.needsUpdate = true;
      updateCap();
      updateStatus();
      view.invalidate();
    }

    // ---- Fade out the p_x ≤ 0 points ("D" or legend button) ---------------------------
    let hollowFaded = false;
    let cancelHollowFade = null;
    function setHollowFaded(value) {
      hollowFaded = value;
      hollowButton?.setAttribute("aria-pressed", String(hollowFaded));
      legend?.classList.toggle("hh-hollow-faded", hollowFaded);
      const from = groups[0]?.hollow.material.uniforms.uOpacity.value ?? 1;
      const to = hollowFaded ? SETTINGS.hollowFadedOpacity : 1;
      cancelHollowFade?.();
      if (hollowButton) hollowButton.textContent = `${hollowFaded ? "Show" : "Hide"} pₓ ≤ 0 (D)`;
      cancelHollowFade = animate(600, eased => {
        for (const group of groups) {
          group.hollow.material.uniforms.uOpacity.value = from + (to - from) * eased;
          group.hollow.visible = group.hollow.material.uniforms.uOpacity.value > 0.01;
        }
        view.invalidate();
      });
      updateStatus();   
    }
    function updateStatus() {
      status.textContent =
        "Filled dots: crossings with pₓ > 0, rings: pₓ ≤ 0. Drag to rotate; scroll or use the slider to move the plane. " +
        `D: ${hollowFaded ? "show" : "hide"} rings · T: ${hint ? "hide" : "show"} trajectories.`;
    }
    setHollowFaded(true);   // initial state is faded px <= 0 
    
    // ---- Controls ---------------------------------------------------------------------
    slider.min = String(xMin);
    slider.max = String(xMax);
    slider.step = step ? String(step) : "any";
    slider.addEventListener("input", () => setPlane(Number(slider.value)));
    canvas.addEventListener("wheel", event => {
      event.preventDefault();
      const delta = event.deltaY * (xMax - xMin) * 0.0002;
      setPlane(c + (step ? Math.sign(delta) * Math.max(step, Math.abs(delta)) : delta));
    }, { passive: false });
    let cancelEntry = null;
    attachRotation(view, canvas, {
      onStart: () => { cancelEntry?.(); cancelEntry = null; lockedSide = null; },
    });
    hollowButton?.addEventListener("click", () => setHollowFaded(!hollowFaded));
    registerSlideKey(canvas, { keyCode: 68, key: "D", description: "Show/hide pₓ ≤ 0 crossings (section slide)" },
      () => setHollowFaded(!hollowFaded));
    registerSlideKey(canvas, { keyCode: 84, key: "T", description: "Show/hide trajectories (section slide)" },
      () => setHint(!hint));

    // ---- Slide entry: continue from slide 2's pose, swing to the default view ---------
    function sideFor(pose) {
      if (SETTINGS.keepSide !== "auto") return SETTINGS.keepSide;
      return Math.sin(pose.azimuth) * Math.sin(pose.elevation) >= 0 ? "below" : "above";
    }

    onSlide(canvas, {
      onEnter: event => {
        view.setActive(true);
        const fromSlide2 = event.previousSlide?.querySelector?.("#hh-3d-plot");
        const target = SETTINGS.defaultCamera;
        if (fromSlide2 || !event.previousSlide) {
          const targetSide = sideFor(target);
          setPlane(SETTINGS.planeStart === "edge"
            ? (targetSide === "below" ? xMax : xMin)
            : SETTINGS.planeStart, { jump: true });
          if (fromSlide2 && handover.pose) {
            const start = { ...handover.pose };
            // Shortest way round for the azimuth.
            const turns = Math.round((start.azimuth - target.azimuth) / (2 * Math.PI));
            const targetAzimuth = target.azimuth + 2 * Math.PI * turns;
            Object.assign(camera, { azimuth: start.azimuth, elevation: start.elevation, zoom: start.zoom });
            lockedSide = targetSide;
            cancelEntry?.();
            cancelEntry = animate(SETTINGS.entryAnimationMs, eased => {
              camera.azimuth = start.azimuth + (targetAzimuth - start.azimuth) * eased;
              camera.elevation = start.elevation + (target.elevation - start.elevation) * eased;
              camera.zoom = start.zoom + (target.zoom - start.zoom) * eased;
              view.invalidate();
            }, () => { cancelEntry = null; lockedSide = null; view.invalidate(); });
          } else {
            Object.assign(camera, target);
          }
        }
        view.invalidate();
      },
      onLeave: () => {
        cancelEntry?.();
        cancelEntry = null;
        lockedSide = null;
        view.setActive(false);
      },
    });

    if (DEBUG) {
      // Debug helpers (only with ?debug in the URL): set the camera / plane from the console.
      window.__hhSetCamera = (azimuth, elevation, zoom = camera.zoom) => {
        Object.assign(camera, { azimuth, elevation, zoom });
        view.invalidate();
      };
      window.__hhSetPlane = value => setPlane(value, { jump: true });
    }

    // Warm-up: the first crossing computation compiles the hot loop (~40 ms); do it
    // during loading so the first plane move is already fast.
    for (const group of groups) crossingsAt(0, group.entry.file, group.buffer);

    // updateStatus();
    setPlane(c);
  } catch (error) {
    console.error(error);
    status.textContent = error.message;
  }
}

// ===================================================================================
// Slide 4 – Poincaré map at x = 0 (former code; data from poincare-map/map/*.json)
// ===================================================================================
const numberedColors = [
  "#f48484", "#7ec8c8", "#e8c97d", "#ae87f7",
  "#9ed07a", "#f0a56b", "#7ea8e8", "#e88ec5",
];
const mapRoot = "../../data/henon-heiles/presentation/poincare-map/";

function readSectionMap(data, filename) {
  if (!Array.isArray(data.y) || !Array.isArray(data.py) ||
      data.y.length !== data.py.length || !data.y.length) {
    throw new Error(`${filename} must contain matching finite y and py arrays.`);
  }
  for (let index = 0; index < data.y.length; index += 1) {
    if (!Number.isFinite(data.y[index]) || !Number.isFinite(data.py[index])) {
      throw new Error(`${filename} contains a non-finite section point.`);
    }
  }
  return data.y.map((y, index) => ({ x: 0, y, py: data.py[index] }));
}

async function loadSectionMaps() {
  const maps = await Promise.all(numberedColors.map(async (color, index) => {
    const filename = `map/${index + 1}.json`;
    const response = await fetch(new URL(`${mapRoot}${filename}`, import.meta.url));
    if (response.status === 404) return null;
    if (!response.ok) throw new Error(`Could not load ${filename} (${response.status}).`);
    return { color, points: readSectionMap(await response.json(), filename) };
  }));
  return maps.filter(Boolean);
}

function getBounds(points) {
  const bounds = { y: [Infinity, -Infinity], py: [Infinity, -Infinity] };
  for (const point of points) {
    for (const coordinate of ["y", "py"]) {
      bounds[coordinate][0] = Math.min(bounds[coordinate][0], point[coordinate]);
      bounds[coordinate][1] = Math.max(bounds[coordinate][1], point[coordinate]);
    }
  }
  return bounds;
}

function prepareCanvas(canvas) {
  const width = canvas.clientWidth;
  const height = canvas.clientHeight;
  if (!width || !height) return null;
  // Sharp at any size: account for reveal.js scaling, not only devicePixelRatio.
  const pixelRatio = Math.min((window.devicePixelRatio || 1) * revealScale(), 3);
  const pixelWidth = Math.round(width * pixelRatio);
  const pixelHeight = Math.round(height * pixelRatio);
  if (canvas.width !== pixelWidth || canvas.height !== pixelHeight) {
    canvas.width = pixelWidth;
    canvas.height = pixelHeight;
  }
  const context = canvas.getContext("2d");
  context.setTransform(pixelRatio, 0, 0, pixelRatio, 0, 0);
  context.clearRect(0, 0, width, height);
  return { context, width, height };
}

async function initializeFixedSection() {
  const fixedCanvas = document.querySelector("#hh-fixed-section-plot");
  const fixedStatus = document.querySelector("#hh-fixed-section-status");
  if (!fixedCanvas || !fixedStatus) return;
  try {
    const maps = await loadSectionMaps();
    if (!maps.length) {
      fixedStatus.textContent =
        "Map data is not present yet. Add numbered JSON files (1.json–8.json) with t, y, and py arrays under data/henon-heiles/presentation/poincare-map/map/.";
      return;
    }
    const bounds = getBounds(maps.flatMap(map => map.points));
    const view = { zoom: 1, panX: 0, panY: 0 };
    let frame = 0;

    function draw() {
      frame = 0;
      const prepared = prepareCanvas(fixedCanvas);
      if (!prepared) return;
      const { context, width, height } = prepared;
      const spanY = Math.max(bounds.y[1] - bounds.y[0], 1e-9);
      const spanPy = Math.max(bounds.py[1] - bounds.py[0], 1e-9);
      const scale = Math.min((width - 72) / spanY, (height - 96) / spanPy) * view.zoom;
      const centerY = (bounds.y[0] + bounds.y[1]) / 2;
      const centerPy = (bounds.py[0] + bounds.py[1]) / 2;
      for (const map of maps) {
        context.beginPath();
        for (const point of map.points) {
          const x = width / 2 + view.panX + (point.y - centerY) * scale;
          const y = height / 2 + view.panY - (point.py - centerPy) * scale;
          context.moveTo(x + 1.8, y);
          context.arc(x, y, 1.8, 0, Math.PI * 2);
        }
        context.fillStyle = map.color;
        context.fill();
      }
      drawSectionMapAxes(context, height - 40);
    }
    const render = () => { if (!frame) frame = requestAnimationFrame(draw); };

    fixedCanvas.addEventListener("wheel", event => {
      event.preventDefault();
      const rect = fixedCanvas.getBoundingClientRect();
      const s = revealScale();
      const pointerX = (event.clientX - rect.left) / s - fixedCanvas.clientWidth / 2;
      const pointerY = (event.clientY - rect.top) / s - fixedCanvas.clientHeight / 2;
      const previousZoom = view.zoom;
      view.zoom = clamp(view.zoom * Math.exp(-event.deltaY * 0.001), 0.5, 24);
      const zoomRatio = view.zoom / previousZoom;
      view.panX = pointerX - (pointerX - view.panX) * zoomRatio;
      view.panY = pointerY - (pointerY - view.panY) * zoomRatio;
      render();
    }, { passive: false });
    let dragPosition;
    fixedCanvas.addEventListener("pointerdown", event => {
      if (event.button !== 0) return;
      dragPosition = { x: event.clientX, y: event.clientY };
      fixedCanvas.setPointerCapture(event.pointerId);
      fixedCanvas.style.cursor = "grabbing";
    });
    fixedCanvas.addEventListener("pointermove", event => {
      if (!dragPosition) return;
      const s = revealScale();
      view.panX += (event.clientX - dragPosition.x) / s;
      view.panY += (event.clientY - dragPosition.y) / s;
      dragPosition = { x: event.clientX, y: event.clientY };
      render();
    });
    const endDrag = () => {
      dragPosition = undefined;
      fixedCanvas.style.cursor = "grab";
    };
    fixedCanvas.addEventListener("pointerup", endDrag);
    fixedCanvas.addEventListener("pointercancel", endDrag);
    fixedCanvas.addEventListener("lostpointercapture", endDrag);
    fixedCanvas.style.cursor = "grab";
    new ResizeObserver(render).observe(fixedCanvas);
    whenReveal().then(reveal => reveal?.on("resize", render));
    render();
    fixedStatus.textContent = "Drag to pan; scroll over the map to zoom.";
  } catch (error) {
    fixedStatus.textContent = error.message;
  }
}

initializeMovingSection();
initializeFixedSection();
