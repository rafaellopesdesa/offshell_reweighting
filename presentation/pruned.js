/* Validated fixed-data restriction to 54 ranked Gaussian auxiliary nuisances. */
(() => {
  'use strict';
  const d = window.PRUNED_RESULTS;
  if (!d || !d.complete || d.points !== 42 || d.profiled_nuisances !== 57) throw new Error('Complete pruned scan required.');
  const fmt = (x, n = 4) => Number(x).toFixed(n);
  const row = (id, values) => {
    const tr = document.createElement('tr');
    for (const value of values) { const td = document.createElement('td'); td.textContent = value; tr.appendChild(td); }
    document.querySelector(id + ' tbody').appendChild(tr);
  };
  const value = (id, text) => { document.querySelector(id).textContent = text; };
  value('#pruned-mu', fmt(d.mu_hat, 8));
  value('#pruned-rms', fmt(d.comparison.rms_delta_q, 5));
  value('#pruned-qzero', fmt(d.q_zero, 5));
  value('#pruned-qtwo', fmt(d.q_two, 5));
  value('#pruned-validation-note', `All 42 profiles converge. Maximum independent ΔNLL disagreement: ${Number(d.checks.maximum_profile_delta_nll_difference).toExponential(2)}. Maximum SR integral mismatch: ${Number(d.checks.maximum_sr_integral_difference).toExponential(2)}.`);
  const names = {superseded_build_validation: 'Superseded setup/check', build_validation: 'Build + construction checks', global_fit: 'Global fit + recovery check', scan_array: '42 independent profiles', finalize: 'Merge + likelihood checks'};
  const duration = x => `${Math.floor(x / 60)}m ${String(Math.round(x % 60)).padStart(2, '0')}s`;
  for (const r of d.resources.phases) row('#pruned-resources', [names[r.phase], r.jobs, duration(r.wall_seconds), fmt(r.cpu_seconds / 60, 2), fmt(r.largest_job_peak_rss_gib, 2) + ' GiB']);
  const c = d.profile_cost;
  value('#pruned-cost-note', `Per-point minimizer wall time: ${fmt(c.fit_seconds.minimum, 1)}–${fmt(c.fit_seconds.maximum, 1)} s (median ${fmt(c.fit_seconds.median, 1)} s). Process peak RSS: ${fmt(c.process_peak_rss_gib.minimum, 2)}–${fmt(c.process_peak_rss_gib.maximum, 2)} GiB, including workspace loading.`);
  value('#snapshot-stamp', 'COMPLETE · ' + d.created_utc.replace('T', ' ').replace('+00:00', ' UTC'));
})();
