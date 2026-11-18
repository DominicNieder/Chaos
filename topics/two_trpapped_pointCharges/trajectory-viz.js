// Slides 1 and 2: configuration space (top view) and 3D view (x, y, p_y) of the
// chaotic and the regular trajectory, with a shared timeline. Slide 2 ends with a
// crossfade from the long chaotic orbit to the solid energy shell (fragment step).

import {
  THREE, animate, attachRotation, clamp, createLineMaterial, createPointMaterial,
  createPolyline, createShellMaterial, createTrajectoryLine, createView, fitPointsFromGeometry,
  handover, insetsFor, loadGroup, loadManifest, onSlide, registerSlideKey, shellGeometry, whenReveal,
} from "./hh-scene.js";
import { zeroVelocityCurve } from "./hh-energy.js";

// ---- Settings (tune on the real projector) -------------------------------------------
const LINE_WIDTH = 1.7;            // CSS px of the slide
const LONG_CHAOTIC_LINE_WIDTH = 0.8;   // the long chaotic orbit fills the shell densely
const LINE_DEPTH_MIN = 0.4;        // brightness of the farthest line parts
const LONG_CHAOTIC_DEPTH_MIN = 0.22;
const MARKER_SIZE = 11;
const CROSSFADE_MS = 1500;
const SHELL_TOGGLE_MS = 300;
const ENTRY_ELEVATION = 0.78;
const ENTRY_ROTATION_MS = 1400;
const ZVC_COLOR = "#9aa3b2";

const configCanvas = document.querySelector("#hh-configuration-plot");
const threeDCanvas = document.querySelector("#hh-3d-plot");
const configStatus = document.querySelector("#hh-viz-status");
const threeDStatus = document.querySelector("#hh-3d-status");
const timelines = [document.querySelector("#hh-config-timeline"), document.querySelector("#hh-3d-timeline")];
const timeLabels = [document.querySelector("#hh-config-time"), document.querySelector("#hh-3d-time")];

function formatTick(value) {
  return Number(value.toPrecision(3)).toString();
}

function setStatus(element, text) {
  if (element) element.textContent = text;
}

async function initialize() {
  if (!configCanvas || !threeDCanvas) return;
  const manifest = await loadManifest();
  const E = manifest.energy;
  const entries = await loadGroup("1trajectory");
  const find = (kind, duration) => entries.find(e => e.kind === kind && e.duration === duration);
  const sets = {
    short: [find("chaotic", "short"), find("regular", "short")].filter(Boolean),
    long: [find("chaotic", "long"), find("regular", "long")].filter(Boolean),
  };
  if (!sets.short.length || !sets.long.length) throw new Error("Short and long trajectories are required.");

  const shell = shellGeometry(E);
  const fitPoints = fitPointsFromGeometry(shell);
  // Zero-velocity curve V(x, y) = E drawn at p_y = 0: it is the boundary of
  // configuration space in the top view and lies exactly on the shell's equator.
  const zvc2d = zeroVelocityCurve(E);
  const zvc = new Float32Array((zvc2d.length / 2) * 3);
  for (let i = 0; i < zvc2d.length / 2; i += 1) zvc.set([zvc2d[2 * i], zvc2d[2 * i + 1], 0], 3 * i);

  let progress = 0;

  // ---- Build one view (slide 1: fixed top view; slide 2: rotatable 3D) ----------------
  function buildView(canvas, camera, { withShell }) {
    const container = canvas.closest(".hh-viz");
    const view = createView({
      container,
      canvas,
      axisCanvas: container.querySelector(".hh-axis-overlay"),
      camera,
      fitPoints,
      insets: insetsFor(container),
    });
    const lines = { short: [], long: [] };
    const markers = { short: [], long: [] };
    for (const duration of ["short", "long"]) {
      for (const entry of sets[duration]) {
        const dense = duration === "long" && entry.kind === "chaotic";
        const material = view.trackLineMaterial(createLineMaterial({
          color: entry.color,
          linewidth: dense ? LONG_CHAOTIC_LINE_WIDTH : LINE_WIDTH,
          depthMin: dense ? LONG_CHAOTIC_DEPTH_MIN : LINE_DEPTH_MIN,
        }));
        const line = view.add(createTrajectoryLine(entry.file, material));
        line.renderOrder = 1;
        lines[duration].push({ entry, line, material });

        const geometry = new THREE.BufferGeometry();
        geometry.setAttribute("position", new THREE.BufferAttribute(new Float32Array(3), 3));
        const pointMaterial = view.trackPointMaterial(createPointMaterial({
          color: entry.color, size: MARKER_SIZE, style: "marker",
        }));
        pointMaterial.depthTest = false;
        const marker = view.add(new THREE.Points(geometry, pointMaterial));
        marker.frustumCulled = false;
        marker.renderOrder = 3;
        markers[duration].push({ entry, marker });
      }
    }
    const zvcMaterial = view.trackLineMaterial(createLineMaterial({ color: ZVC_COLOR, linewidth: 1.2, depthMin: 0.6 }));
    view.add(createPolyline(zvc, zvcMaterial));

    let shellMesh = null;
    if (withShell) {
      shellMesh = view.add(new THREE.Mesh(shell, createShellMaterial()));
      shellMesh.visible = false;
      shellMesh.renderOrder = 0;
    }
    return { view, lines, markers, shellMesh };
  }

  const configCamera = { azimuth: 0, elevation: 0, zoom: 1 };
  const threeDCamera = { azimuth: 0, elevation: 0, zoom: 1 };
  const config = buildView(configCanvas, configCamera, { withShell: false });
  const threeD = buildView(threeDCanvas, threeDCamera, { withShell: true });

  // ---- Timeline -------------------------------------------------------------------------
  function applyProgress(target) {
    const duration = progress >= 1 ? "long" : "short";
    for (const key of ["short", "long"]) {
      target.lines[key].forEach(({ entry, line }, index) => {
        const visible = key === duration;
        line.visible = visible && lineVisibility(target);
        const count = entry.file.count;
        const last = key === "long" || progress >= 1 ? count - 1 : Math.floor(clamp(progress, 0, 1) * (count - 1));
        line.setDrawCount(last + 1);
        const { marker } = target.markers[key][index];
        marker.visible = visible && line.visible;
        const p = entry.file.positions;
        marker.geometry.attributes.position.array.set(p.subarray(3 * last, 3 * last + 3));
        marker.geometry.attributes.position.needsUpdate = true;
      });
    }
    target.view.invalidate();
  }

  function currentTimeLabel() {
    if (progress >= 1) return `t = ${formatTick(sets.long[0].file.meta.tEnd)} a.u.`;
    return `t = ${formatTick(progress * sets.short[0].file.meta.tEnd)} a.u.`;
  }

  function setProgress(value) {
    progress = clamp(value, 0, 1);
    timelines.forEach(timeline => {
      if (timeline) timeline.value = String(Math.round(progress * Number(timeline.max)));
    });
    timeLabels.forEach(label => { if (label) label.textContent = currentTimeLabel(); });
    applyProgress(config);
    applyProgress(threeD);
  }

  timelines.forEach(timeline => timeline?.addEventListener("input", () => {
    setProgress(Number(timeline.value) / Number(timeline.max));
    const message = progress >= 1 ? "Long-term trajectory" :
      `Short trajectory · t = ${formatTick(progress * sets.short[0].file.meta.tEnd)}`;
    setStatus(configStatus, message);
    setStatus(threeDStatus, message);
  }));

  // ---- Slide 2: shell state machine ------------------------------------------------------
  // solid: set by the fragment step (crossfade lines -> solid shell)
  // shellOn: toggled with "E"; without `solid` the shell is a cutaway (back faces only)
  const shellState = { solid: false, shellOn: false, shellOnBeforeSolid: false };
  let fade = 0;                     // 0 = lines, 1 = solid shell (crossfade value)
  let cancelFade = null;
  const shellMaterial = threeD.shellMesh.material;

  function lineVisibility(target) {
    return target !== threeD || fade < 1 || !shellState.shellOn;
  }

  function applyShell() {
    const { solid, shellOn } = shellState;
    const mesh = threeD.shellMesh;
    if (!shellOn) {
      mesh.visible = fade > 0 && fade < 1;
    } else {
      mesh.visible = true;
    }
    const crossfading = fade > 0 && fade < 1;
    if (crossfading) {
      shellMaterial.side = THREE.FrontSide;
      shellMaterial.transparent = true;
      shellMaterial.depthWrite = false;
      shellMaterial.opacity = fade;
    } else {
      shellMaterial.transparent = false;
      shellMaterial.depthWrite = true;
      shellMaterial.opacity = 1;
      shellMaterial.side = solid ? THREE.FrontSide : THREE.BackSide;
    }
    shellMaterial.needsUpdate = true;
    for (const key of ["short", "long"]) {
      for (const { material } of threeD.lines[key]) {
        const lineOpacity = shellOn ? 1 - fade : 1;
        material.transparent = lineOpacity < 1;
        material.opacity = lineOpacity;
        material.needsUpdate = true;
      }
    }
    applyProgress(threeD);
  }

  function runFade(target, duration) {
    cancelFade?.();
    const from = fade;
    cancelFade = animate(duration, eased => {
      fade = from + (target - from) * eased;
      applyShell();
    }, () => { fade = target; cancelFade = null; applyShell(); });
  }

  function setSolid(solid, { immediate = false } = {}) {
    if (solid === shellState.solid) return;
    shellState.solid = solid;
    if (solid) {
      shellState.shellOnBeforeSolid = shellState.shellOn;
      shellState.shellOn = true;
      if (progress < 1) setProgress(1);
    } else {
      shellState.shellOn = shellState.shellOnBeforeSolid;
    }
    if (immediate) { fade = solid ? 1 : 0; applyShell(); } else runFade(solid ? 1 : 0, CROSSFADE_MS);
    updateThreeDStatus();
  }

  function toggleShell() {
    shellState.shellOn = !shellState.shellOn;
    if (shellState.solid) {
      // In the solid end state "E" hides the shell and shows the orbits again (and back).
      runFade(shellState.shellOn ? 1 : 0, SHELL_TOGGLE_MS);
    } else {
      applyShell();
    }
    updateThreeDStatus();
  }

  function updateThreeDStatus(prefix = "") {
    const shellHint = shellState.solid
      ? (shellState.shellOn ? "Energy shell: the orbits fill this region. E: show orbits." : "E: show energy shell.")
      : (shellState.shellOn ? "E: hide energy shell." : "E: show energy shell.");
    setStatus(threeDStatus, `${prefix}Drag to rotate; scroll to zoom; scrub the timeline. ${shellHint}`);
  }

  const fragment = document.querySelector(".hh-shell-fragment");
  whenReveal().then(reveal => {
    if (!reveal || !fragment) return;
    reveal.on("fragmentshown", event => { if (event.fragment === fragment) setSolid(true); });
    reveal.on("fragmenthidden", event => { if (event.fragment === fragment) setSolid(false); });
  });

  registerSlideKey(threeDCanvas, { keyCode: 69, key: "E", description: "Energy shell on/off (Hénon–Heiles slides)" }, toggleShell);

  // ---- Interaction --------------------------------------------------------------------
  let cancelRotation = null;
  attachRotation(threeD.view, threeDCanvas, { onStart: () => { cancelRotation?.(); cancelRotation = null; } });
  threeDCanvas.addEventListener("wheel", event => {
    event.preventDefault();
    threeDCamera.zoom = clamp(threeDCamera.zoom * Math.exp(-event.deltaY * 0.001), 0.65, 2.5);
    threeD.view.invalidate();
  }, { passive: false });

  onSlide(configCanvas, {
    onEnter: () => {
      cancelRotation?.();
      threeDCamera.elevation = 0;
      setProgress(0);
      config.view.setActive(true);
      setStatus(configStatus, "Use the timeline to advance the trajectory.");
    },
    onLeave: () => config.view.setActive(false),
  });

  onSlide(threeDCanvas, {
    onEnter: event => {
      threeD.view.setActive(true);
      setSolid(Boolean(fragment?.classList.contains("visible")), { immediate: true });
      const fromSlide3 = event.previousSlide && event.previousSlide.querySelector?.("#hh-moving-section-plot");
      if (!fromSlide3) {
        setProgress(1);
        const start = threeDCamera.elevation;
        cancelRotation?.();
        cancelRotation = animate(ENTRY_ROTATION_MS, eased => {
          threeDCamera.elevation = start + (ENTRY_ELEVATION - start) * eased;
          threeD.view.invalidate();
        }, () => { cancelRotation = null; });
      }
      updateThreeDStatus();
    },
    onLeave: () => {
      cancelRotation?.();
      cancelRotation = null;
      handover.pose = { ...threeDCamera, solid: shellState.solid && shellState.shellOn };
      threeD.view.setActive(false);
    },
  });

  setProgress(progress);
  applyShell();
}

initialize().catch(error => {
  console.error(error);
  setStatus(configStatus, error.message);
  setStatus(threeDStatus, error.message);
});
