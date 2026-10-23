/* Current six-model family; earlier snapshots remain available on later slides. */
(() => {
  'use strict';
  const data = window.FAMILY_RESULTS;
  if (!data?.complete || data.models.length !== 6 || data.grid.length !== 61) {
    throw new Error('Six complete, validated 61-point scans are required.');
  }
  const fmt = (x, digits) => Number(x).toFixed(digits);
  for (const model of data.models) {
    if (model.points !== 61) throw new Error('Incomplete workspace scan.');
    const section = document.querySelector(`[data-family-model="${model.label}"]`);
    section.querySelector('[data-family-mu]').textContent = fmt(model.mu_hat, 7);
    section.querySelector('[data-family-rms]').textContent = fmt(model.rms_delta_q, 5);
    section.querySelector('[data-family-max]').textContent = fmt(model.max_abs_delta_q, 5);
    section.querySelector('[data-family-max-mu]').textContent = fmt(model.mu_at_max_abs_delta_q, 2);
  }
  document.querySelector('#snapshot-stamp').textContent = '6 SCANS VALIDATED · ' +
    new Date(data.created_utc).toISOString().slice(0, 16).replace('T', ' ') + ' UTC';
  const reading = document.body.classList.contains('reading');
  document.querySelectorAll('[data-family-start]').forEach(link => {
    link.href = reading ? '#family-true-4096' : '#/family-true-4096';
  });
  if (reading) {
    const scrollToFamily = () => {
      if (!location.hash.startsWith('#/family-')) return;
      const id = location.hash.slice(2);
      const section = document.getElementById(id);
      if (section) { history.replaceState(null, '', '#' + id); section.scrollIntoView(); }
    };
    window.addEventListener('load', scrollToFamily);
    window.addEventListener('hashchange', scrollToFamily);
  }
})();
