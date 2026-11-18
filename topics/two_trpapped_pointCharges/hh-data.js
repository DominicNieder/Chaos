import { buildXIndex } from "./hh-energy.js";

const manifestUrl = new URL(
  "../../data/henon-heiles/presentation/bin/manifest.json",
  import.meta.url,
);
const trajectoryCache = new Map();
let manifestPromise;

async function fetchResource(url, description) {
  const response = await fetch(url);
  if (!response.ok) throw new Error(`Could not load ${description} (${response.status}).`);
  return response.arrayBuffer();
}

function decodeFloat32(buffer, count) {
  if (buffer.byteLength !== count * 4 * Float32Array.BYTES_PER_ELEMENT) {
    throw new Error("Trajectory coordinate file has an invalid byte length.");
  }
  const view = new DataView(buffer);
  const values = new Float32Array(count * 4);
  for (let index = 0; index < values.length; index += 1) {
    values[index] = view.getFloat32(index * 4, true);
  }
  return values;
}

function decodeFloat64(buffer, count) {
  if (buffer.byteLength !== count * Float64Array.BYTES_PER_ELEMENT) {
    throw new Error("Trajectory time file has an invalid byte length.");
  }
  const view = new DataView(buffer);
  const values = new Float64Array(count);
  for (let index = 0; index < count; index += 1) {
    values[index] = view.getFloat64(index * 8, true);
  }
  return values;
}

function validateTrajectoryEntry(entry) {
  if (!entry || typeof entry.id !== "string" ||
      !Number.isInteger(entry.count) || entry.count < 2 ||
      typeof entry.coordinatesFile !== "string") {
    throw new Error("The trajectory manifest contains an invalid entry.");
  }
}

export async function loadManifest() {
  manifestPromise ??= (async () => {
    const response = await fetch(manifestUrl);
    if (!response.ok) throw new Error(`Could not load trajectory manifest (${response.status}).`);
    const manifest = await response.json();
    if (manifest.format !== "henon-heiles-trajectories" || manifest.version !== 1 ||
        manifest.coordinateType !== "float32-le" || manifest.timeType !== "float64-le" ||
        !Number.isFinite(manifest.energy) || !Array.isArray(manifest.trajectories)) {
      throw new Error("Trajectory manifest format or version is unsupported.");
    }
    for (const entry of manifest.trajectories) validateTrajectoryEntry(entry);
    return manifest;
  })();
  return manifestPromise;
}

async function loadTrajectory(entry) {
  const key = `${entry.coordinatesFile}:${entry.timesFile ?? `${entry.t0}:${entry.dt}`}`;
  if (!trajectoryCache.has(key)) {
    trajectoryCache.set(key, (async () => {
      const coordinatesUrl = new URL(entry.coordinatesFile, manifestUrl);
      const coordinateBuffer = await fetchResource(coordinatesUrl, entry.coordinatesFile);
      const samples = decodeFloat32(coordinateBuffer, entry.count);
      let times;
      if (entry.timesFile) {
        const timesUrl = new URL(entry.timesFile, manifestUrl);
        const timeBuffer = await fetchResource(timesUrl, entry.timesFile);
        times = decodeFloat64(timeBuffer, entry.count);
      } else {
        if (!Number.isFinite(entry.t0) || !Number.isFinite(entry.dt) || entry.dt <= 0) {
          throw new Error(`Trajectory ${entry.id} has invalid uniform time metadata.`);
        }
        times = Float64Array.from(
          { length: entry.count },
          (_, index) => entry.t0 + index * entry.dt,
        );
      }

      for (let index = 0; index < times.length; index += 1) {
        if (!Number.isFinite(times[index]) ||
            (index > 0 && times[index] <= times[index - 1])) {
          throw new Error(`Trajectory ${entry.id} has non-increasing sample times.`);
        }
      }
      return { samples, times, xIndex: undefined };
    })());
  }
  return { ...entry, ...(await trajectoryCache.get(key)) };
}

export async function loadTrajectories(ids, { withXIndex = false } = {}) {
  const manifest = await loadManifest();
  const entries = ids.map(id => {
    const entry = manifest.trajectories.find(candidate => candidate.id === id);
    if (!entry) throw new Error(`Trajectory "${id}" is missing from the manifest.`);
    return entry;
  });
  const trajectories = await Promise.all(entries.map(loadTrajectory));
  if (withXIndex) {
    for (const trajectory of trajectories) {
      trajectory.xIndex ??= buildXIndex(trajectory.samples);
    }
  }
  return { manifest, trajectories };
}
