export const AXIS_COLORS = Object.freeze({
  x: "#ff7b72",
  y: "#7ee787",
  py: "#79c0ff",
});

function drawAxisArrows(context, height, axes) {
  const origin = { x: 28, y: height - 58 };
  for (const axis of axes) {
    const length = Math.hypot(axis.direction.x, axis.direction.y);
    if (length < 0.08) continue;

    const dx = axis.direction.x / length * 36;
    const dy = -axis.direction.y / length * 36;
    const endX = origin.x + dx;
    const endY = origin.y + dy;
    const angle = Math.atan2(dy, dx);
    const headLength = 8;

    context.beginPath();
    context.moveTo(origin.x, origin.y);
    context.lineTo(endX, endY);
    context.strokeStyle = axis.color;
    context.lineWidth = 3.5;
    context.lineCap = "round";
    context.stroke();

    context.beginPath();
    context.moveTo(endX, endY);
    context.lineTo(
      endX - headLength * Math.cos(angle - Math.PI / 6),
      endY - headLength * Math.sin(angle - Math.PI / 6),
    );
    context.lineTo(
      endX - headLength * Math.cos(angle + Math.PI / 6),
      endY - headLength * Math.sin(angle + Math.PI / 6),
    );
    context.closePath();
    context.fillStyle = axis.color;
    context.fill();
    context.font = "bold 14px sans-serif";
    context.fillText(axis.label, endX + 5, endY - 4);
  }
}

export function drawCameraAxes(context, height, camera) {
  const cosAzimuth = Math.cos(camera.azimuth);
  const sinAzimuth = Math.sin(camera.azimuth);
  const cosElevation = Math.cos(camera.elevation);
  const sinElevation = Math.sin(camera.elevation);
  drawAxisArrows(context, height, [
    { direction: { x: cosAzimuth, y: sinAzimuth * cosElevation }, label: "x", color: AXIS_COLORS.x },
    { direction: { x: -sinAzimuth, y: cosAzimuth * cosElevation }, label: "y", color: AXIS_COLORS.y },
    { direction: { x: 0, y: -sinElevation }, label: "pᵧ", color: AXIS_COLORS.py },
  ]);
}

export function drawSectionMapAxes(context, height) {
  drawAxisArrows(context, height, [
    { direction: { x: 1, y: 0 }, label: "y", color: AXIS_COLORS.y },
    { direction: { x: 0, y: 1 }, label: "pᵧ", color: AXIS_COLORS.py },
  ]);
}
