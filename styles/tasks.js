document.addEventListener("DOMContentLoaded", () => {
  const key = "tasks-" + window.location.pathname;
  const boxes = document.querySelectorAll(".task-list-item-checkbox");
  const saved = JSON.parse(localStorage.getItem(key) || "{}");
  boxes.forEach((cb, i) => {
    cb.disabled = false;
    if (saved[i]) cb.checked = true;
    cb.addEventListener("change", () => {
      const state = {};
      boxes.forEach((b, j) => (state[j] = b.checked));
      localStorage.setItem(key, JSON.stringify(state));
    });
  });
});
