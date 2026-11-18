const ESCAPE_ENERGY = 1 / 6;

export function V(x, y) {
  return 0.5 * (x * x + y * y) + x * x * y - (y * y * y) / 3;
}

export function H(state) {
  const [x, y, px, py] = state;
  return 0.5 * (px * px + py * py) + V(x, y);
}

export function f(state) {
  const [x, y, px, py] = state;
  return [px, py, -x - 2 * x * y, -y - x * x + y * y];
}

function assertBoundedEnergy(energy) {
  if (!Number.isFinite(energy) || energy <= 0 || energy >= ESCAPE_ENERGY) {
    throw new RangeError(`energy must be in (0, ${ESCAPE_ENERGY})`);
  }
}

function radialBoundary(theta, energy) {
  const sinTheta = Math.sin(theta);
  const cubic = sinTheta - (4 / 3) * sinTheta ** 3;
  let high;
  if (cubic < 0) {
    high = -1 / (3 * cubic);
    if (0.5 * high ** 2 + cubic * high ** 3 <= energy) {
      throw new Error("No bounded radial boundary exists at this angle.");
    }
  } else {
    high = Math.max(1, Math.sqrt(2 * energy));
    while (0.5 * high ** 2 + cubic * high ** 3 < energy) high *= 2;
  }

  let low = 0;
  for (let iteration = 0; iteration < 64; iteration += 1) {
    const middle = (low + high) / 2;
    const value = 0.5 * middle ** 2 + cubic * middle ** 3;
    if (value < energy) low = middle;
    else high = middle;
  }
  return (low + high) / 2;
}

export function zeroVelocityCurve(energy, resolution = 256) {
  assertBoundedEnergy(energy);
  const count = Math.max(16, Math.floor(resolution));
  const points = new Float32Array((count + 1) * 2);
  for (let index = 0; index <= count; index += 1) {
    const theta = (2 * Math.PI * index) / count;
    const radius = radialBoundary(theta, energy);
    points[index * 2] = radius * Math.cos(theta);
    points[index * 2 + 1] = radius * Math.sin(theta);
  }
  return points;
}

function shellResolution(resolution) {
  if (typeof resolution === "number") {
    const angular = Math.max(16, Math.floor(resolution));
    return { angular, radial: Math.max(8, Math.floor(angular / 3)) };
  }
  const angular = Math.max(16, Math.floor(resolution?.angular ?? 144));
  const radial = Math.max(8, Math.floor(resolution?.radial ?? 48));
  return { angular, radial };
}

export function buildShellMesh(energy, resolution = {}) {
  assertBoundedEnergy(energy);
  const { angular, radial } = shellResolution(resolution);
  const verticesPerSheet = 1 + angular * radial;
  const positions = new Float32Array(verticesPerSheet * 2 * 3);
  const normals = new Float32Array(positions.length);
  const indices = [];
  let minX = Infinity;
  let maxX = -Infinity;
  let minY = Infinity;
  let maxY = -Infinity;
  let maxPy = Math.sqrt(2 * energy);

  for (const [sheet, sign] of [[0, 1], [1, -1]]) {
    const vertexOffset = sheet * verticesPerSheet;
    const centerPy = sign * Math.sqrt(2 * energy);
    positions.set([0, 0, centerPy], vertexOffset * 3);
    normals.set([0, 0, sign], vertexOffset * 3);

    for (let ring = 1; ring <= radial; ring += 1) {
      const fraction = ring / radial;
      for (let angleIndex = 0; angleIndex < angular; angleIndex += 1) {
        const theta = (2 * Math.PI * angleIndex) / angular;
        const radius = fraction * radialBoundary(theta, energy);
        const x = radius * Math.cos(theta);
        const y = radius * Math.sin(theta);
        const py = sign * Math.sqrt(Math.max(0, 2 * (energy - V(x, y))));
        const vertex = vertexOffset + 1 + (ring - 1) * angular + angleIndex;
        const positionOffset = vertex * 3;
        positions.set([x, y, py], positionOffset);

        const normal = [x + 2 * x * y, y + x * x - y * y, py];
        const normalLength = Math.hypot(...normal);
        normals.set(normal.map(value => value / normalLength), positionOffset);
        minX = Math.min(minX, x);
        maxX = Math.max(maxX, x);
        minY = Math.min(minY, y);
        maxY = Math.max(maxY, y);
        maxPy = Math.max(maxPy, Math.abs(py));

        const indexInRing = angleIndex;
        const nextAngle = (angleIndex + 1) % angular;
        const indexOf = (ringIndex, angle) =>
          vertexOffset + 1 + (ringIndex - 1) * angular + angle;
        if (ring === 1) {
          if (sign > 0) {
            indices.push(vertexOffset, indexOf(1, indexInRing), indexOf(1, nextAngle));
          } else {
            indices.push(vertexOffset, indexOf(1, nextAngle), indexOf(1, indexInRing));
          }
        } else {
          const innerCurrent = indexOf(ring - 1, indexInRing);
          const innerNext = indexOf(ring - 1, nextAngle);
          const outerCurrent = indexOf(ring, indexInRing);
          const outerNext = indexOf(ring, nextAngle);
          if (sign > 0) {
            indices.push(
              innerCurrent, outerCurrent, outerNext,
              innerCurrent, outerNext, innerNext,
            );
          } else {
            indices.push(
              innerCurrent, outerNext, outerCurrent,
              innerCurrent, innerNext, outerNext,
            );
          }
        }
      }
    }
  }

  return {
    positions,
    normals,
    indices: new Uint32Array(indices),
    bounds: { x: [minX, maxX], y: [minY, maxY], py: [-maxPy, maxPy] },
  };
}

export function shellBounds(energy) {
  const bounds = buildShellMesh(energy).bounds;
  const extent = shellXExtent(energy);
  return { ...bounds, x: [-extent, extent] };
}

export function shellXExtent(energy) {
  let low = 0;
  let high = 1;
  while (capInterval(high, energy)) high *= 2;
  for (let iteration = 0; iteration < 64; iteration += 1) {
    const middle = (low + high) / 2;
    if (capInterval(middle, energy)) low = middle;
    else high = middle;
  }
  return (low + high) / 2;
}

function bisectRoot(fn, low, high, lowValue, highValue) {
  if (lowValue === 0) return low;
  if (highValue === 0) return high;
  if (Math.sign(lowValue) === Math.sign(highValue)) {
    throw new Error("Root bracket does not straddle zero.");
  }
  for (let iteration = 0; iteration < 80; iteration += 1) {
    const middle = (low + high) / 2;
    const middleValue = fn(middle);
    if (Math.abs(middleValue) < 1e-14 || high - low < 1e-13) return middle;
    if (Math.sign(middleValue) === Math.sign(lowValue)) {
      low = middle;
      lowValue = middleValue;
    } else {
      high = middle;
      highValue = middleValue;
    }
  }
  return (low + high) / 2;
}

export function capInterval(x, energy) {
  assertBoundedEnergy(energy);
  if (!Number.isFinite(x)) throw new TypeError("x must be finite.");

  const root = Math.sqrt(1 + 4 * x * x);
  const localMinimum = (1 - root) / 2;
  const localMaximum = (1 + root) / 2;
  const value = y => V(x, y) - energy;
  const minimumValue = value(localMinimum);
  const maximumValue = value(localMaximum);
  if (minimumValue > 1e-13 || maximumValue <= 0) return null;
  if (Math.abs(minimumValue) <= 1e-13) return [localMinimum, localMinimum];

  let left = localMinimum - 1;
  let step = 1;
  while (value(left) < 0 && step < 1e6) {
    step *= 2;
    left = localMinimum - step;
  }
  if (value(left) < 0) return null;

  const lower = bisectRoot(value, left, localMinimum, value(left), minimumValue);
  const upper = bisectRoot(value, localMinimum, localMaximum, minimumValue, maximumValue);
  return [lower, upper];
}

export function buildCapGeometry(x, energy, resolution = 128) {
  const interval = capInterval(x, energy);
  if (!interval) {
    return {
      interval: null,
      positions: new Float32Array(),
      indices: new Uint32Array(),
      outline: new Float32Array(),
    };
  }

  const [low, high] = interval;
  if (high - low < 1e-12) {
    return {
      interval,
      positions: new Float32Array([x, low, 0]),
      indices: new Uint32Array(),
      outline: new Float32Array([x, low, 0]),
    };
  }

  const count = Math.max(8, Math.floor(resolution));
  const positions = new Float32Array((count + 1) * 2 * 3);
  const outline = new Float32Array((count + 1) * 2 * 3);
  const indices = [];
  const yValues = [];
  const pyValues = [];
  for (let index = 0; index <= count; index += 1) {
    const u = (1 - Math.cos(Math.PI * index / count)) / 2;
    const y = low + u * (high - low);
    const py = Math.sqrt(Math.max(0, 2 * (energy - V(x, y))));
    yValues.push(y);
    pyValues.push(py);

    const vertex = index * 2;
    positions.set([x, y, -py, x, y, py], vertex * 3);
    outline.set([x, y, -py], index * 3);
    outline.set([x, y, py], (2 * count + 1 - index) * 3);

    if (index < count) {
      const next = vertex + 2;
      indices.push(vertex, next, next + 1, vertex, next + 1, vertex + 1);
    }
  }

  return {
    interval,
    yValues: new Float64Array(yValues),
    pyValues: new Float64Array(pyValues),
    positions,
    indices: new Uint32Array(indices),
    outline,
  };
}

export function buildXIndex(samples, bins = 256) {
  if (!(samples instanceof Float32Array || samples instanceof Float64Array) ||
      samples.length < 8 || samples.length % 4 !== 0) {
    throw new TypeError("samples must be an interleaved Float32Array or Float64Array of states.");
  }
  const count = samples.length / 4;
  const binCount = Math.max(1, Math.floor(bins));
  let minX = Infinity;
  let maxX = -Infinity;
  for (let index = 0; index < count; index += 1) {
    const x = samples[index * 4];
    minX = Math.min(minX, x);
    maxX = Math.max(maxX, x);
  }
  const buckets = Array.from({ length: binCount }, () => []);
  const binFor = x => Math.min(
    binCount - 1,
    Math.max(0, Math.floor(((x - minX) / (maxX - minX || 1)) * binCount)),
  );
  for (let index = 0; index < count - 1; index += 1) {
    const x0 = samples[index * 4];
    const x1 = samples[(index + 1) * 4];
    const firstBin = binFor(Math.min(x0, x1));
    const lastBin = binFor(Math.max(x0, x1));
    for (let bin = firstBin; bin <= lastBin; bin += 1) buckets[bin].push(index);
  }
  return { minX, maxX, bins: buckets, binCount };
}

function hermiteCoefficients(samples, index, duration) {
  const offset = index * 4;
  const nextOffset = offset + 4;
  const start = samples.subarray(offset, nextOffset);
  const end = samples.subarray(nextOffset, nextOffset + 4);
  const startDerivative = f(start);
  const endDerivative = f(end);
  return Array.from({ length: 4 }, (_, coordinate) => {
    const first = start[coordinate];
    const second = end[coordinate];
    const firstSlope = duration * startDerivative[coordinate];
    const secondSlope = duration * endDerivative[coordinate];
    return [
      2 * first - 2 * second + firstSlope + secondSlope,
      -3 * first + 3 * second - 2 * firstSlope - secondSlope,
      firstSlope,
      first,
    ];
  });
}

function evaluateCubic(coefficients, fraction) {
  return ((coefficients[0] * fraction + coefficients[1]) * fraction +
    coefficients[2]) * fraction + coefficients[3];
}

function derivativeCubic(coefficients, fraction) {
  return (3 * coefficients[0] * fraction + 2 * coefficients[1]) * fraction +
    coefficients[2];
}

function segmentCrosses(x0, x1, x) {
  if (x1 > x0) return x0 <= x && x < x1;
  if (x1 < x0) return x1 <= x && x < x0;
  return false;
}

export function crossingsAt(x, samples, times, index, output = null) {
  if (!Number.isFinite(x)) throw new TypeError("section x must be finite.");
  if (!(times instanceof Float64Array) || times.length * 4 !== samples.length) {
    throw new TypeError("times must be a Float64Array matching the sample count.");
  }
  if (output !== null && !(output instanceof Float32Array)) {
    throw new TypeError("output must be a Float32Array.");
  }
  if (x < index.minX || x > index.maxX) return output ? 0 : new Float32Array();

  const bin = Math.min(
    index.binCount - 1,
    Math.max(0, Math.floor(((x - index.minX) / (index.maxX - index.minX || 1)) * index.binCount)),
  );
  const candidates = index.bins[bin];
  const crossings = output ? null : [];
  let count = 0;
  let previousTime = -Infinity;
  for (const sampleIndex of candidates) {
    const offset = sampleIndex * 4;
    const x0 = samples[offset];
    const x1 = samples[offset + 4];
    if (!segmentCrosses(x0, x1, x)) continue;
    const duration = times[sampleIndex + 1] - times[sampleIndex];
    if (!(duration > 0)) throw new Error(`Non-increasing time at sample ${sampleIndex}.`);

    const coefficients = hermiteCoefficients(samples, sampleIndex, duration);
    const xCoefficients = coefficients[0];
    let low = 0;
    let high = 1;
    let fraction = (x - x0) / (x1 - x0);
    const increasing = x1 > x0;
    const xAt = value => evaluateCubic(xCoefficients, value);
    for (let iteration = 0; iteration < 2; iteration += 1) {
      const residual = xAt(fraction) - x;
      if (Math.abs(residual) < 1e-13) break;
      if (increasing ? residual < 0 : residual > 0) low = fraction;
      else high = fraction;
      const derivative = derivativeCubic(xCoefficients, fraction);
      const candidate = fraction - residual / derivative;
      fraction = Number.isFinite(candidate) && candidate > low && candidate < high
        ? candidate
        : (low + high) / 2;
    }

    for (let iteration = 0; iteration < 32; iteration += 1) {
      const residual = xAt(fraction) - x;
      if (Math.abs(residual) < 1e-13 || high - low < 1e-13) break;
      if (increasing ? residual < 0 : residual > 0) low = fraction;
      else high = fraction;
      fraction = (low + high) / 2;
    }

    const state = coefficients.map(coordinate => evaluateCubic(coordinate, fraction));
    const crossingTime = times[sampleIndex] + fraction * duration;
    if (Math.abs(crossingTime - previousTime) <=
        1e-12 * Math.max(1, Math.abs(crossingTime))) continue;
    if (output) {
      if (3 * count + 2 >= output.length) {
        throw new RangeError("output buffer is too small for the section crossings.");
      }
      output[3 * count] = state[1];
      output[3 * count + 1] = state[3];
      output[3 * count + 2] = state[2];
    } else {
      crossings.push(state[1], state[3], state[2]);
    }
    count += 1;
    previousTime = crossingTime;
  }
  return output ? count : new Float32Array(crossings);
}
