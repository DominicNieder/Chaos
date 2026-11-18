import assert from "node:assert/strict";
import { readFile } from "node:fs/promises";
import { dirname, join, resolve } from "node:path";
import { fileURLToPath } from "node:url";
import test from "node:test";

import {
  H,
  V,
  buildCapGeometry,
  buildShellMesh,
  buildXIndex,
  capInterval,
  crossingsAt,
  zeroVelocityCurve,
} from "../../topics/two_trpapped_pointCharges/hh-energy.js";

const root = resolve(dirname(fileURLToPath(import.meta.url)), "../..");
const dataRoot = join(root, "data/henon-heiles/presentation");
const binaryRoot = join(dataRoot, "bin");
const manifest = JSON.parse(await readFile(join(binaryRoot, "manifest.json"), "utf8"));
const energy = manifest.energy;

async function loadTrajectory(entry) {
  const coordinateBuffer = await readFile(join(binaryRoot, entry.coordinatesFile));
  const samples = new Float32Array(entry.count * 4);
  const coordinateView = new DataView(
    coordinateBuffer.buffer,
    coordinateBuffer.byteOffset,
    coordinateBuffer.byteLength,
  );
  for (let index = 0; index < samples.length; index += 1) {
    samples[index] = coordinateView.getFloat32(index * 4, true);
  }

  let times;
  if (entry.timesFile) {
    const timeBuffer = await readFile(join(binaryRoot, entry.timesFile));
    const timeView = new DataView(timeBuffer.buffer, timeBuffer.byteOffset, timeBuffer.byteLength);
    times = new Float64Array(entry.count);
    for (let index = 0; index < times.length; index += 1) {
      times[index] = timeView.getFloat64(index * 8, true);
    }
  } else {
    times = Float64Array.from(
      { length: entry.count },
      (_, index) => entry.t0 + index * entry.dt,
    );
  }
  return { samples, times };
}

function trajectoryEntry(id) {
  const entry = manifest.trajectories.find(trajectory => trajectory.id === id);
  assert.ok(entry, `manifest contains ${id}`);
  return entry;
}

test("positive crossings at x=0 agree with the saved chaotic map", async () => {
  const entry = trajectoryEntry("trajectory-1");
  const { samples, times } = await loadTrajectory(entry);
  const index = buildXIndex(samples);
  const crossings = crossingsAt(0, samples, times, index);
  const capacity = Math.max(...index.bins.map(bucket => bucket.length));
  const output = new Float32Array(3 * capacity);
  const outputCount = crossingsAt(0, samples, times, index, output);
  assert.equal(outputCount, crossings.length / 3);
  assert.deepEqual(output.slice(0, 3 * outputCount), crossings);
  const positive = [];
  for (let index = 0; index < crossings.length; index += 3) {
    if (crossings[index + 2] > 0) positive.push([crossings[index], crossings[index + 1]]);
  }
  const map = JSON.parse(await readFile(join(dataRoot, "poincare-map/map/1.json"), "utf8"));

  assert.equal(positive.length, map.y.length);
  for (let index = 0; index < positive.length; index += 1) {
    assert.ok(Math.abs(positive[index][0] - map.y[index]) <= 1e-4);
    assert.ok(Math.abs(positive[index][1] - map.py[index]) <= 1e-4);
  }
});

test("crossings from every trajectory lie inside each sampled cap", async () => {
  const entries = manifest.trajectories.filter(entry => entry.id.startsWith("trajectory-"));
  const cValues = [-0.2, 0, 0.2];
  for (const entry of entries) {
    const { samples, times } = await loadTrajectory(entry);
    const index = buildXIndex(samples);
    for (const c of cValues) {
      const interval = capInterval(c, energy);
      assert.ok(interval, `cap exists at x=${c}`);
      const crossings = crossingsAt(c, samples, times, index);
      for (let offset = 0; offset < crossings.length; offset += 3) {
        const [y, py] = [crossings[offset], crossings[offset + 1]];
        assert.ok(y >= interval[0] - 1e-5 && y <= interval[1] + 1e-5);
        assert.ok(0.5 * py * py + V(c, y) <= energy + 1e-6);
      }
    }
  }
});

test("shell and cap boundaries satisfy the energy equation", () => {
  const shell = buildShellMesh(energy, { angular: 96, radial: 32 });
  let maxResidual = 0;
  for (let index = 0; index < shell.positions.length; index += 3) {
    const residual = H([
      shell.positions[index],
      shell.positions[index + 1],
      0,
      shell.positions[index + 2],
    ]) - energy;
    maxResidual = Math.max(maxResidual, Math.abs(residual));
  }
  assert.ok(maxResidual <= 5e-5, `shell residual ${maxResidual}`);

  for (const c of [-0.2, 0, 0.2]) {
    const cap = buildCapGeometry(c, energy, 128);
    assert.ok(cap.interval);
    for (let index = 0; index < cap.outline.length; index += 3) {
      const [x, y, py] = cap.outline.slice(index, index + 3);
      assert.ok(Math.abs(x - c) <= 1e-6);
      assert.ok(Math.abs(0.5 * py * py + V(x, y) - energy) <= 1e-6);
    }
  }

  const curve = zeroVelocityCurve(energy, 256);
  for (let index = 0; index < curve.length; index += 2) {
    assert.ok(Math.abs(V(curve[index], curve[index + 1]) - energy) <= 5e-5);
  }
});
