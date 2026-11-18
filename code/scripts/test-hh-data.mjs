import assert from "node:assert/strict";
import { readFile } from "node:fs/promises";
import { dirname, join, resolve } from "node:path";
import { fileURLToPath } from "node:url";
import test from "node:test";

const root = resolve(dirname(fileURLToPath(import.meta.url)), "../..");
const originalFetch = globalThis.fetch;
globalThis.fetch = async input => {
  const url = new URL(input);
  const bytes = await readFile(fileURLToPath(url));
  return new Response(bytes);
};

const { loadTrajectories } = await import(
  "../../topics/two_trpapped_pointCharges/hh-data.js"
);

test("binary loader supports long and replacement trajectories with distinct metadata", async () => {
  const { trajectories: long } = await loadTrajectories([
    "chaotic-long",
    "regular-long",
  ]);
  const [chaotic, regular] = long;

  assert.equal(chaotic.count, 214914);
  assert.equal(regular.count, 10825);
  assert.equal(chaotic.samples.length, chaotic.count * 4);
  assert.equal(regular.samples.length, regular.count * 4);
  assert.ok(chaotic.times.at(-1) > regular.times.at(-1));
  assert.ok(chaotic.times[1] > chaotic.times[0]);
  assert.ok(regular.times[1] > regular.times[0]);

  const { trajectories: sourceMap } = await loadTrajectories(["trajectory-1"]);
  assert.equal(sourceMap[0].id, "trajectory-1");
  assert.equal(sourceMap[0].color, "#f48484");
  assert.equal(sourceMap[0].label, "Chaotic");
  assert.strictEqual(sourceMap[0].samples, chaotic.samples);
  assert.strictEqual(sourceMap[0].times, chaotic.times);
});

test.after(() => {
  globalThis.fetch = originalFetch;
});
