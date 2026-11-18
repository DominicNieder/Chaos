// Shared WebGL (three.js) infrastructure for the Hénon–Heiles slides.
//
// World coordinates: X = x, Y = y, Z = p_y.
// Camera model (same as the former Canvas2D code), orthographic:
//   e1 = ( cos a,       −sin a,        0     )   screen right
//   e2 = ( sin a cos e,  cos a cos e,  −sin e )   screen up
//   e3 = ( sin a sin e,  cos a sin e,   cos e )   towards the viewer (e3 = e1 × e2)
// screenX = x cos a − y sin a, screenY = r_y cos e − p_y sin e, depth = r_y sin e + p_y cos e.

import * as THREE from "../../vendor/three/three.module.js";
import { Line2 } from "../../vendor/three/lines/Line2.js";
import { LineGeometry } from "../../vendor/three/lines/LineGeometry.js";
import { LineMaterial } from "../../vendor/three/lines/LineMaterial.js";
import { drawCameraAxes } from "./axis-indicator.js";
import { loadManifest as loadDataManifest, loadTrajectories } from "./hh-data.js";
import { buildShellMesh, crossingsAt as interpolateCrossings } from "./hh-energy.js";

export { THREE };

export const DEBUG = new URLSearchParams(window.location.search).has("debug");
const MAX_PIXEL_RATIO = 3;

// ---- Data loading (cached per page) ------------------------------------------------
const cache = new Map();

function cached(key, factory) {
  if (!cache.has(key)) cache.set(key, factory());
  return cache.get(key);
}

export const loadManifest = () => cached("manifest", loadDataManifest);

function renderFile(trajectory) {
  const { samples, times, xIndex } = trajectory;
  const count = times.length;
  const positions = new Float32Array(3 * count);
  for (let index = 0; index < count; index += 1) {
    positions[3 * index] = samples[4 * index];
    positions[3 * index + 1] = samples[4 * index + 1];
    positions[3 * index + 2] = samples[4 * index + 3];
  }
  let geometry = null;
  let indexSummary = null;
  return {
    key: trajectory.id,
    meta: { ...trajectory.bounds, tEnd: times[count - 1] },
    state: samples,
    times,
    count,
    positions,
    trajectory,
    get index() {
      if (!xIndex) return null;
      indexSummary ??= {
        ...xIndex,
        maxBin: Math.max(1, ...xIndex.bins.map(bucket => bucket.length)),
      };
      return indexSummary;
    },
    lineGeometry() {
      if (!geometry) {
        geometry = new LineGeometry();
        geometry.setPositions(positions);
      }
      return geometry;
    },
  };
}

/** Adapt original manifest IDs to the renderer's visual groups. */
export async function loadGroup(name) {
  const idsByGroup = {
    "1trajectory": ["chaotic-short", "regular-short", "chaotic-long", "regular-long"],
    "poincare-map": Array.from({ length: 8 }, (_, index) => `trajectory-${index + 1}`),
  };
  const ids = idsByGroup[name];
  if (!ids) throw new Error(`Unknown trajectory group "${name}".`);
  const { trajectories } = await loadTrajectories(ids, { withXIndex: name === "poincare-map" });
  return trajectories.map(trajectory => ({
    ...trajectory,
    kind: trajectory.id.startsWith("chaotic") || trajectory.id === "trajectory-1"
      ? "chaotic"
      : "regular",
    duration: trajectory.id.endsWith("-short") ? "short" : "long",
    file: renderFile(trajectory),
  }));
}

export function crossingsAt(x, file, output) {
  const trajectory = file.trajectory;
  return interpolateCrossings(x, trajectory.samples, trajectory.times, trajectory.xIndex, output);
}

/** Shell mesh geometry for energy E (shared, built once). */
export function shellGeometry(E) {
  return cached(`shell:${E}`, () => {
    const mesh = buildShellMesh(E);
    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute("position", new THREE.BufferAttribute(mesh.positions, 3));
    geometry.setAttribute("normal", new THREE.BufferAttribute(mesh.normals, 3));
    geometry.setIndex(new THREE.BufferAttribute(mesh.indices, 1));
    geometry.computeBoundingSphere();
    return geometry;
  });
}

/** A subsample of shell vertices used to fit the view tightly for every camera pose. */
export function fitPointsFromGeometry(geometry, maxPoints = 600) {
  const source = geometry.getAttribute("position").array;
  const n = source.length / 3;
  const step = Math.max(1, Math.floor(n / maxPoints));
  const points = [];
  for (let i = 0; i < n; i += step) points.push(source[3 * i], source[3 * i + 1], source[3 * i + 2]);
  return new Float32Array(points);
}

// ---- Camera ------------------------------------------------------------------------
export function cameraBasis(azimuth, elevation) {
  const ca = Math.cos(azimuth);
  const sa = Math.sin(azimuth);
  const ce = Math.cos(elevation);
  const se = Math.sin(elevation);
  return {
    e1: new THREE.Vector3(ca, -sa, 0),
    e2: new THREE.Vector3(sa * ce, ca * ce, -se),
    e3: new THREE.Vector3(sa * se, ca * se, ce),
  };
}

export function clamp(value, min, max) {
  return Math.max(min, Math.min(max, value));
}

export function smoothstep(t) {
  return t * t * (3 - 2 * t);
}

// ---- Materials ---------------------------------------------------------------------
/**
 * Fat line material with depth cueing: brightness falls from 1 (nearest) to
 * depthMin (farthest) across the depth range set by the view.
 */
export function createLineMaterial({ color, linewidth = 1.6, depthMin = 0.35, clippingPlanes = null }) {
  const material = new LineMaterial({ color, linewidth, worldUnits: false });
  material.clippingPlanes = clippingPlanes;
  material.userData.depthUniforms = {
    uHhNear: { value: 0 },
    uHhFar: { value: 1 },
    uHhDepthMin: { value: depthMin },
  };
  material.onBeforeCompile = shader => {
    Object.assign(shader.uniforms, material.userData.depthUniforms);
    shader.vertexShader = "varying float vHhDepth;\n" + shader.vertexShader.replace(
      "#include <clipping_planes_vertex>",
      "#include <clipping_planes_vertex>\n\tvHhDepth = -mvPosition.z;",
    );
    shader.fragmentShader =
      "varying float vHhDepth;\nuniform float uHhNear;\nuniform float uHhFar;\nuniform float uHhDepthMin;\n" +
      shader.fragmentShader.replace(
        "#include <tonemapping_fragment>",
        "float hhDepth = clamp((vHhDepth - uHhNear) / max(uHhFar - uHhNear, 1e-6), 0.0, 1.0);\n" +
        "\tgl_FragColor.rgb *= mix(1.0, uHhDepthMin, hhDepth);\n\t#include <tonemapping_fragment>",
      );
  };
  material.customProgramCacheKey = () => "hh-depth-cued-line";
  return material;
}

/** Round point sprites: filled discs, hollow rings, or discs with a white rim (marker). */
export function createPointMaterial({ color, size = 4, style = "filled", brightness = 1, offset = 0 }) {
  return new THREE.ShaderMaterial({
    transparent: true,
    depthWrite: false,
    uniforms: {
      uColor: { value: new THREE.Color(color) },
      uSize: { value: size },
      uPixelRatio: { value: 1 },
      uBrightness: { value: brightness },
      uOpacity: { value: 1 },
      uStyle: { value: style === "hollow" ? 1 : style === "marker" ? 2 : 0 },
      uToward: { value: new THREE.Vector3(0, 0, 1) },
      uOffset: { value: offset },
    },
    vertexShader: /* glsl */`
      uniform float uSize;
      uniform float uPixelRatio;
      uniform vec3 uToward;
      uniform float uOffset;
      void main() {
        vec3 p = position + uToward * uOffset;
        gl_Position = projectionMatrix * modelViewMatrix * vec4(p, 1.0);
        gl_PointSize = uSize * uPixelRatio;
      }
    `,
    fragmentShader: /* glsl */`
      uniform vec3 uColor;
      uniform float uBrightness;
      uniform float uOpacity;
      uniform float uStyle;
      uniform float uSize;
      uniform float uPixelRatio;
      void main() {
        vec2 q = gl_PointCoord * 2.0 - 1.0;
        float r = length(q);
        float aa = 2.0 / max(uSize * uPixelRatio, 1.0);
        float alpha = 1.0 - smoothstep(1.0 - aa, 1.0, r);
        vec3 color = uColor * uBrightness;
        if (uStyle > 0.5 && uStyle < 1.5) {
          alpha *= smoothstep(0.52 - aa, 0.52 + aa, r);   // ring
        } else if (uStyle > 1.5) {
          color = mix(color, vec3(1.0), smoothstep(0.62 - aa, 0.62 + aa, r));   // white rim
        }
        alpha *= uOpacity;
        if (alpha < 0.02) discard;
        gl_FragColor = vec4(color, alpha);
        #include <colorspace_fragment>
      }
    `,
  });
}

export function createShellMaterial({ color = "#3d4b5e", side = THREE.FrontSide, clippingPlanes = null } = {}) {
  return new THREE.MeshStandardMaterial({
    color,
    roughness: 0.62,
    metalness: 0.05,
    side,
    clippingPlanes,
  });
}

/** Line2 drawn from a shared LineGeometry; drawCount = number of samples to show. */
export function createTrajectoryLine(file, material) {
  const line = new Line2(file.lineGeometry(), material);
  line.frustumCulled = false;
  line.setDrawCount = samples => {
    // Line2 is instanced (one instance per segment). Shared geometry, so the
    // count is applied per draw call via onBeforeRender.
    line.userData.segments = Math.max(0, samples - 1);
  };
  line.userData.segments = file.count - 1;
  line.onBeforeRender = (renderer, scene, camera, geometry) => {
    geometry.instanceCount = line.userData.segments;
  };
  return line;
}

/** A Line2 for a small polyline that changes (cap outline, zero-velocity curve). */
export function createPolyline(positions, material) {
  const geometry = new LineGeometry();
  geometry.setPositions(positions);
  const line = new Line2(geometry, material);
  line.frustumCulled = false;
  line.setPositions = next => {
    line.geometry.dispose();
    line.geometry = new LineGeometry();
    line.geometry.setPositions(next);
  };
  return line;
}

// ---- reveal.js integration -----------------------------------------------------------
let revealPromise = null;

/** Resolves with the Reveal API once it is ready, or with null outside reveal.js. */
export function whenReveal() {
  if (!revealPromise) {
    revealPromise = new Promise(resolve => {
      const attempt = () => {
        const reveal = window.Reveal;
        if (!reveal) return false;
        if (reveal.isReady?.()) resolve(reveal);
        else reveal.on("ready", () => resolve(reveal));
        return true;
      };
      if (attempt()) return;
      const onLoad = () => { if (!attempt()) resolve(null); };
      if (document.readyState === "complete") setTimeout(onLoad, 0);
      else window.addEventListener("load", onLoad, { once: true });
    });
  }
  return revealPromise;
}

export function revealScale() {
  const scale = window.Reveal?.getScale?.();
  return Number.isFinite(scale) && scale > 0 ? scale : 1;
}

/**
 * Call onEnter(event) / onLeave(event) when the slide containing `element`
 * becomes current / stops being current. Outside reveal.js, onEnter is called once.
 */
export async function onSlide(element, { onEnter, onLeave }) {
  const reveal = await whenReveal();
  const slideOf = node => node?.closest?.("section") ?? null;
  const mySlide = slideOf(element);
  if (!reveal) { onEnter?.({}); return; }
  let current = false;
  const update = event => {
    const isCurrent = event.currentSlide === mySlide;
    if (isCurrent && !current) onEnter?.(event);
    if (!isCurrent && current) onLeave?.(event);
    current = isCurrent;
  };
  update({ currentSlide: reveal.getCurrentSlide(), previousSlide: null });
  reveal.on("slidechanged", update);
}

// Custom keys: one reveal.js binding per key, dispatched to the handler of the
// current slide (several slides may use the same key).
const keyHandlers = new Map();   // keyCode -> [{ slide, handler, description }]

export async function registerSlideKey(element, { keyCode, key, description }, handler) {
  const reveal = await whenReveal();
  const slide = element.closest("section");
  if (!keyHandlers.has(keyCode)) {
    keyHandlers.set(keyCode, []);
    const dispatch = () => {
      const current = reveal ? reveal.getCurrentSlide() : slide;
      for (const entry of keyHandlers.get(keyCode)) {
        if (entry.slide === current) entry.handler();
      }
    };
    if (reveal?.addKeyBinding) {
      reveal.addKeyBinding({ keyCode, key, description }, dispatch);
    } else {
      window.addEventListener("keydown", event => {
        if (event.keyCode === keyCode && !event.ctrlKey && !event.metaKey && !event.altKey) dispatch();
      });
    }
  }
  keyHandlers.get(keyCode).push({ slide, handler, description });
}

// Camera hand-over between slides (slide 2 end pose -> slide 3 start pose).
export const handover = { pose: null };

// ---- View ------------------------------------------------------------------------------
/**
 * One WebGL view in a slide container.
 * options: { container, canvas, axisCanvas, camera: {azimuth, elevation, zoom},
 *            fitPoints (Float32Array xyz), insets: () => ({top, right, bottom, left}),
 *            padding }
 */
export function createView({ container, canvas, axisCanvas = null, camera, fitPoints, insets, padding = 1.06 }) {
  const renderer = new THREE.WebGLRenderer({ canvas, antialias: true, alpha: true });
  renderer.setClearColor(0x000000, 0);
  renderer.localClippingEnabled = true;
  const scene = new THREE.Scene();
  const ortho = new THREE.OrthographicCamera(-1, 1, 1, -1, 0.01, 100);
  scene.add(ortho);

  // Lights follow the camera, so the shell is lit the same way from every view.
  scene.add(new THREE.HemisphereLight(0xdfe8ff, 0x1a1d24, 1.1));
  const keyLight = new THREE.DirectionalLight(0xffffff, 1.6);
  keyLight.position.set(-0.6, 0.8, 1.0);
  ortho.add(keyLight);
  ortho.add(keyLight.target);
  keyLight.target.position.set(0, 0, -1);

  const lineMaterials = new Set();
  const pointMaterials = new Set();
  const listeners = { beforeRender: [] };
  let active = false;
  let frameRequested = false;
  let size = { width: 0, height: 0 };
  let pixelRatio = 1;
  const depthRange = { near: 0, far: 1 };
  const DISTANCE = 10;

  function resize() {
    const width = container.clientWidth;
    const height = container.clientHeight;
    if (!width || !height) return false;
    const ratio = Math.min((window.devicePixelRatio || 1) * revealScale(), MAX_PIXEL_RATIO);
    if (width !== size.width || height !== size.height || ratio !== pixelRatio) {
      size = { width, height };
      pixelRatio = ratio;
      renderer.setPixelRatio(ratio);
      renderer.setSize(width, height, false);
      if (axisCanvas) {
        const w = axisCanvas.clientWidth;
        const h = axisCanvas.clientHeight;
        axisCanvas.width = Math.round(w * ratio);
        axisCanvas.height = Math.round(h * ratio);
      }
    }
    return true;
  }

  function updateCamera() {
    const { e1, e2, e3 } = cameraBasis(camera.azimuth, camera.elevation);
    // Tight fit: project the fit points onto the screen axes for this pose.
    let minX = Infinity; let maxX = -Infinity;
    let minY = Infinity; let maxY = -Infinity;
    let minD = Infinity; let maxD = -Infinity;
    for (let i = 0; i < fitPoints.length; i += 3) {
      const x = fitPoints[i]; const y = fitPoints[i + 1]; const z = fitPoints[i + 2];
      const sx = x * e1.x + y * e1.y + z * e1.z;
      const sy = x * e2.x + y * e2.y + z * e2.z;
      const sd = x * e3.x + y * e3.y + z * e3.z;
      if (sx < minX) minX = sx; if (sx > maxX) maxX = sx;
      if (sy < minY) minY = sy; if (sy > maxY) maxY = sy;
      if (sd < minD) minD = sd; if (sd > maxD) maxD = sd;
    }
    const inset = insets ? insets() : { top: 0, right: 0, bottom: 0, left: 0 };
    const availableWidth = Math.max(size.width - inset.left - inset.right, 10);
    const availableHeight = Math.max(size.height - inset.top - inset.bottom, 10);
    const scale = Math.min(
      availableWidth / Math.max((maxX - minX) * padding, 1e-9),
      availableHeight / Math.max((maxY - minY) * padding, 1e-9),
    ) * camera.zoom;                                    // CSS px per world unit
    const halfWidth = size.width / (2 * scale);
    const halfHeight = size.height / (2 * scale);
    // Content centre is placed at the centre of the available area.
    const offsetX = (inset.left - inset.right) / 2 / scale;
    const offsetY = (inset.bottom - inset.top) / 2 / scale;
    const cx = (minX + maxX) / 2 - offsetX;
    const cy = (minY + maxY) / 2 - offsetY;
    const centre = new THREE.Vector3()
      .addScaledVector(e1, cx)
      .addScaledVector(e2, cy)
      .addScaledVector(e3, (minD + maxD) / 2);
    ortho.left = -halfWidth; ortho.right = halfWidth;
    ortho.top = halfHeight; ortho.bottom = -halfHeight;
    ortho.near = 0.01; ortho.far = 2 * DISTANCE;
    ortho.position.copy(centre).addScaledVector(e3, DISTANCE);
    ortho.up.copy(e2);
    ortho.lookAt(centre);
    ortho.updateProjectionMatrix();
    ortho.updateMatrixWorld();
    // View-space depth (distance from the camera plane) of the content.
    const centreDepth = (minD + maxD) / 2;
    depthRange.near = DISTANCE - (maxD - centreDepth);
    depthRange.far = DISTANCE + (centreDepth - minD);
    return { e1, e2, e3, scale };
  }

  function render() {
    frameRequested = false;
    if (!active || !resize()) return;
    const start = DEBUG ? performance.now() : 0;
    const basis = updateCamera();
    for (const material of lineMaterials) {
      material.resolution.set(size.width, size.height);
      const u = material.userData.depthUniforms;
      u.uHhNear.value = depthRange.near;
      u.uHhFar.value = depthRange.far;
    }
    for (const material of pointMaterials) {
      material.uniforms.uPixelRatio.value = pixelRatio;
      material.uniforms.uToward.value.copy(basis.e3);
    }
    for (const listener of listeners.beforeRender) listener(basis);
    renderer.render(scene, ortho);
    if (axisCanvas) {
      const context = axisCanvas.getContext("2d");
      context.setTransform(1, 0, 0, 1, 0, 0);
      context.clearRect(0, 0, axisCanvas.width, axisCanvas.height);
      context.setTransform(pixelRatio, 0, 0, pixelRatio, 0, 0);
      context.translate(18, 0);
      drawCameraAxes(context, axisCanvas.clientHeight, camera);
    }
    if (DEBUG) {
      const elapsed = performance.now() - start;
      console.debug(`[hh] ${canvas.id}: frame ${elapsed.toFixed(2)} ms`);
    }
  }

  function invalidate() {
    if (!active || frameRequested) return;
    frameRequested = true;
    requestAnimationFrame(render);
  }

  const resizeObserver = new ResizeObserver(invalidate);
  resizeObserver.observe(container);
  whenReveal().then(reveal => reveal?.on("resize", invalidate));

  return {
    renderer,
    scene,
    camera,
    container,
    get pixelRatio() { return pixelRatio; },
    add(object) { scene.add(object); return object; },
    trackLineMaterial(material) { lineMaterials.add(material); return material; },
    trackPointMaterial(material) { pointMaterials.add(material); return material; },
    onBeforeRender(listener) { listeners.beforeRender.push(listener); },
    setActive(value) {
      active = value;
      if (active) invalidate();
    },
    get active() { return active; },
    invalidate,
    basis: () => cameraBasis(camera.azimuth, camera.elevation),
  };
}

/** Drag-to-rotate with the same feel as the former Canvas2D code. */
export function attachRotation(view, canvas, { minElevation = -0.3, maxElevation = 1.3, onStart } = {}) {
  let drag = null;
  canvas.addEventListener("pointerdown", event => {
    if (event.button !== 0) return;
    drag = { x: event.clientX, y: event.clientY };
    canvas.setPointerCapture(event.pointerId);
    onStart?.();
  });
  canvas.addEventListener("pointermove", event => {
    if (!drag) return;
    view.camera.azimuth += (event.clientX - drag.x) * 0.008;
    view.camera.elevation = clamp(view.camera.elevation + (event.clientY - drag.y) * 0.008, minElevation, maxElevation);
    drag = { x: event.clientX, y: event.clientY };
    view.invalidate();
  });
  const end = () => { drag = null; };
  canvas.addEventListener("pointerup", end);
  canvas.addEventListener("pointercancel", end);
  canvas.addEventListener("lostpointercapture", end);
}

/** Small tween helper driven by requestAnimationFrame. Returns a cancel function. */
export function animate(duration, step, done) {
  const start = performance.now();
  let frame = 0;
  const tick = now => {
    const t = clamp((now - start) / duration, 0, 1);
    step(smoothstep(t), t);
    if (t < 1) frame = requestAnimationFrame(tick);
    else { frame = 0; done?.(); }
  };
  frame = requestAnimationFrame(tick);
  return () => { if (frame) cancelAnimationFrame(frame); frame = 0; };
}

/** Bottom inset for the view fit, from the controls overlay height (CSS px). */
export function insetsFor(container, { top = 16, side = 24 } = {}) {
  return () => {
    const controls = container.querySelector(".hh-viz-controls");
    return { top, right: side, left: side, bottom: (controls?.offsetHeight ?? 0) + 12 };
  };
}
