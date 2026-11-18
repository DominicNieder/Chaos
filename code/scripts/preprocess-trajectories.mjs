import { createHash } from "node:crypto";
import { mkdir, readFile, writeFile } from "node:fs/promises";
import { dirname, join, relative, resolve } from "node:path";
import { fileURLToPath } from "node:url";

const root = resolve(dirname(fileURLToPath(import.meta.url)), "../..");
const dataRoot = join(root, "data/henon-heiles/presentation");
const outputRoot = join(dataRoot, "bin");
const energyTarget = 0.1342;
const parameter = { a: 1, m: 1, w: 1 };
const trajectoryInputs = [
  {
    id: "chaotic-short",
    label: "Chaotic",
    color: "#f48484",
    file: "1trajectory/chaotic_short.json",
    referencedBy: ["slide-1", "slide-2-short"],
  },
  {
    id: "regular-short",
    label: "Regular",
    color: "#69d2e7",
    file: "1trajectory/regular_short.json",
    referencedBy: ["slide-1", "slide-2-short"],
  },
  {
    id: "chaotic-long",
    label: "Chaotic",
    color: "#f48484",
    file: "1trajectory/chaotic_long.json",
    referencedBy: ["slide-2-long"],
  },
  {
    id: "regular-long",
    label: "Regular",
    color: "#69d2e7",
    file: "1trajectory/regular_long.json",
    referencedBy: ["slide-2-long"],
  },
  ...Array.from({ length: 8 }, (_, index) => ({
    id: `trajectory-${index + 1}`,
    label: index === 0 ? "Chaotic" : `Regular trajectory ${index + 1}`,
    color: [
      "#f48484", "#7ec8c8", "#e8c97d", "#ae87f7",
      "#9ed07a", "#f0a56b", "#7ea8e8", "#e88ec5",
    ][index],
    file: `poincare-map/traj/${index + 1}.json`,
    referencedBy: ["slide-3", "slide-4-source"],
  })),
];

function hamiltonian(state) {
  const [x, y, px, py] = state;
  return (px * px + py * py) / (2 * parameter.m) +
    0.5 * parameter.m * parameter.w ** 2 * (x * x + y * y) +
    parameter.a * (x * x * y - (y * y * y) / 3);
}

function isUniform(times) {
  if (times.length < 3) return { uniform: true, dt: 0 };
  const dt = times[1] - times[0];
  const tolerance = 1e-12 * Math.max(1, Math.abs(dt));
  return {
    uniform: times.slice(2).every((time, index) =>
      Math.abs((time - times[index + 1]) - dt) <= tolerance),
    dt,
  };
}

function serializeFloat32LittleEndian(points) {
  const buffer = Buffer.allocUnsafe(points.length * 4 * Float32Array.BYTES_PER_ELEMENT);
  let offset = 0;
  for (const point of points) {
    for (const value of point) {
      buffer.writeFloatLE(value, offset);
      offset += Float32Array.BYTES_PER_ELEMENT;
    }
  }
  return buffer;
}

function serializeFloat64LittleEndian(times) {
  const buffer = Buffer.allocUnsafe(times.length * Float64Array.BYTES_PER_ELEMENT);
  times.forEach((time, index) =>
    buffer.writeDoubleLE(time, index * Float64Array.BYTES_PER_ELEMENT));
  return buffer;
}

function validateSectionMap(data, filename) {
  if (!Array.isArray(data.y) || !Array.isArray(data.py) ||
      data.y.length !== data.py.length ||
      (data.t && data.t.length !== data.y.length)) {
    throw new Error(`${filename} has inconsistent y, py, or t arrays.`);
  }
  if (!data.y.length || !data.y.every(Number.isFinite) ||
      !data.py.every(Number.isFinite) ||
      (data.t && !data.t.every(Number.isFinite))) {
    throw new Error(`${filename} has empty or non-finite section data.`);
  }
  return data.y.length;
}

function formatMiB(bytes) {
  return `${(bytes / (1024 * 1024)).toFixed(2)} MiB`;
}

async function preprocessTrajectory(input, sourceHashes, uniqueFiles) {
  const sourcePath = join(dataRoot, input.file);
  const sourceBytes = await readFile(sourcePath);
  const sourceHash = createHash("sha256").update(sourceBytes).digest("hex");
  const data = JSON.parse(sourceBytes.toString("utf8"));
  if (!Array.isArray(data.t) || !Array.isArray(data.phi) ||
      data.t.length !== data.phi.length || !data.phi.length) {
    throw new Error(`${input.file} must contain matching, non-empty t and phi arrays.`);
  }

  const bounds = {
    x: [Infinity, -Infinity],
    y: [Infinity, -Infinity],
    px: [Infinity, -Infinity],
    py: [Infinity, -Infinity],
  };
  const points = [];
  const times = [];
  let maxEnergyError = 0;
  for (let index = 0; index < data.phi.length; index += 1) {
    const state = data.phi[index];
    const time = data.t[index];
    if (!Array.isArray(state) || state.length < 4 ||
        !Number.isFinite(time) || !state.slice(0, 4).every(Number.isFinite)) {
      throw new Error(`${input.file} contains an invalid state at index ${index}.`);
    }

    maxEnergyError = Math.max(maxEnergyError, Math.abs(hamiltonian(state) - energyTarget));
    for (const [coordinateIndex, coordinate] of ["x", "y", "px", "py"].entries()) {
      bounds[coordinate][0] = Math.min(bounds[coordinate][0], state[coordinateIndex]);
      bounds[coordinate][1] = Math.max(bounds[coordinate][1], state[coordinateIndex]);
    }

    const previousState = data.phi[index - 1];
    const duplicate = index > 0 && time === data.t[index - 1] &&
      state.slice(0, 4).every((value, coordinate) => value === previousState[coordinate]);
    if (!duplicate) {
      points.push(state.slice(0, 4));
      times.push(time);
    }
  }

  if (maxEnergyError > 1e-5) {
    console.warn(`WARNING ${input.file}: max |H - ${energyTarget}| = ${maxEnergyError}`);
  }

  let coordinatesFile;
  let timesFile = null;
  let timeInfo;
  let assetBytes = 0;
  const existing = sourceHashes.get(sourceHash);
  if (existing) {
    coordinatesFile = existing.coordinatesFile;
    timesFile = existing.timesFile;
    timeInfo = existing.timeInfo;
    assetBytes = existing.assetBytes;
  } else {
    const coordinateBytes = serializeFloat32LittleEndian(points);
    const { uniform, dt } = isUniform(times);
    coordinatesFile = `trajectories/${sourceHash}.f32le`;
    await mkdir(dirname(join(outputRoot, coordinatesFile)), { recursive: true });
    await writeFile(join(outputRoot, coordinatesFile), coordinateBytes);
    assetBytes += coordinateBytes.byteLength;

    if (uniform) {
      timeInfo = { t0: times[0], dt };
    } else {
      timesFile = `times/${sourceHash}.f64le`;
      await mkdir(dirname(join(outputRoot, timesFile)), { recursive: true });
      const timeBytes = serializeFloat64LittleEndian(times);
      await writeFile(join(outputRoot, timesFile), timeBytes);
      assetBytes += timeBytes.byteLength;
      timeInfo = { timesFile };
    }
    sourceHashes.set(sourceHash, { coordinatesFile, timesFile, timeInfo, assetBytes });
    uniqueFiles.add(coordinatesFile);
    if (timesFile) uniqueFiles.add(timesFile);
  }

  return {
    ...input,
    source: input.file,
    sourceHash,
    sourceBytes: sourceBytes.byteLength,
    sourcePoints: data.phi.length,
    count: points.length,
    coordinatesFile,
    timesFile,
    ...timeInfo,
    bounds,
    energy: energyTarget,
    maxEnergyError,
    tEnd: times.at(-1),
    outputBytes: assetBytes,
  };
}

async function main() {
  await mkdir(outputRoot, { recursive: true });
  const sourceHashes = new Map();
  const uniqueFiles = new Set();
  const trajectories = [];
  for (const input of trajectoryInputs) {
    trajectories.push(await preprocessTrajectory(input, sourceHashes, uniqueFiles));
  }

  const sectionMaps = [];
  for (let index = 1; index <= 8; index += 1) {
    const file = `poincare-map/map/${index}.json`;
    const data = JSON.parse(await readFile(join(dataRoot, file), "utf8"));
    sectionMaps.push({
      id: `map-${index}`,
      source: file,
      count: validateSectionMap(data, file),
    });
  }

  const manifest = {
    format: "henon-heiles-trajectories",
    version: 1,
    energy: energyTarget,
    coordinateType: "float32-le",
    timeType: "float64-le",
    trajectories: trajectories.map(({ sourceBytes, sourcePoints, outputBytes, ...entry }) => entry),
    sectionMaps,
  };
  await writeFile(
    join(outputRoot, "manifest.json"),
    `${JSON.stringify(manifest, null, 2)}\n`,
  );

  const rows = trajectories.map(entry => [
    entry.source,
    String(entry.sourcePoints),
    String(entry.count),
    formatMiB(entry.sourceBytes),
    formatMiB(entry.outputBytes),
  ]);
  const widths = rows[0].map((_, column) =>
    Math.max(...rows.map(row => row[column].length), ["source", "before", "after", "JSON", "binary"][column].length));
  const header = ["source", "before", "after", "JSON", "binary"];
  console.log([
    header.map((cell, index) => cell.padEnd(widths[index])).join("  "),
    ...rows.map(row => row.map((cell, index) => cell.padEnd(widths[index])).join("  ")),
  ].join("\n"));
  console.log(`Wrote ${uniqueFiles.size} unique binary files and manifest.json under ${relative(root, outputRoot)}.`);
  console.log(`Validated ${sectionMaps.length} map files.`);
}

main().catch(error => {
  console.error(error);
  process.exitCode = 1;
});
