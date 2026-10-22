/* Snapshot data are local and frozen; nothing here queries the cluster. */
(() => {
  'use strict';
  const data = window.RESULTS;
  const reading = window.matchMedia('(max-width:700px)').matches;
  document.body.classList.toggle('reading', reading);
  if (reading) document.querySelector('.masthead > a').href = '#overview';
  const time = data.created_utc.replace('T', ' ').replace('+00:00', ' UTC');
  document.querySelector('#snapshot-stamp').textContent = `${data.status.toUpperCase()} · ${time}`;
  document.querySelectorAll('.snapshot-time').forEach(el => { el.textContent = time; });
  const addRow = (target, values) => {
    const tr = document.createElement('tr');
    values.forEach(value => { const td = document.createElement('td'); td.textContent = value; tr.appendChild(td); });
    document.querySelector(target).appendChild(tr);
  };
  const fmt = (x, places = 4) => Number(x).toFixed(places);
  for (const [size, sample] of Object.entries(data.samples)) {
    document.querySelectorAll(`[data-mu="${size}"]`).forEach(el => { el.textContent = fmt(sample.mu_hat); });
    document.querySelectorAll(`[data-count="${size}"]`).forEach(el => { el.textContent = `${sample.points} / 31 scan points`; });
    document.querySelectorAll(`[data-status="${size}"]`).forEach(el => {
      el.textContent = sample.complete ? 'COMPLETE + VALIDATED' : 'PRELIMINARY';
      el.classList.toggle('complete', sample.complete);
    });
    document.querySelector(`[data-caption="${size}"]`).textContent = sample.complete
      ? 'Completed profile, normalized to the continuous global minimum.'
      : `Partial scan. Markers show converged coordinates only; gaps have not been interpolated. Baseline: ${sample.fit_complete ? 'final global fit' : 'saved best-fit checkpoint'}.`;
    const choices = [.05, .1, .5, 1, 1.5];
    for (const mu of choices) {
      const row = sample.rows.find(r => Math.abs(r.mu - mu) < 1e-9);
      if (row) addRow(`[data-points="${size}"] tbody`, [fmt(mu, 2), fmt(row.q, 3), fmt(row.reference_q, 3)]);
    }
    addRow('#fit-table tbody', [Number(size).toLocaleString('en-US'), ...sample.start_mu_estimates.map(x => fmt(x, 6)),
      sample.start_nll_spread.toExponential(2), sample.fit_complete ? `Quality ${sample.covariance_quality}` : 'Pending']);
    document.querySelectorAll(`[data-shift="${size}"]`).forEach(el => { el.textContent = `+${fmt(100*(sample.mu_hat-1), 3)}%`; });
    const diagnostic = data.quadrature_first_start[size].fitted;
    const corrected = sample.global_diagnostics;
    addRow('#quadrature-table tbody', [Number(size).toLocaleString('en-US'),
      `${fmt(100*diagnostic.sr_integral_relative_difference, 3)}%`,
      `${fmt(100*corrected.sr_integral_relative_difference_at_minimum, 6)}%`,
      fmt(corrected.maximum_absolute_alpha, 5)]);
    addRow('#before-table tbody', [Number(size).toLocaleString('en-US'),
      fmt(sample.before.mu_hat, 6)+(size === '4096' ? '*' : ''), fmt(sample.mu_hat, 6),
      '+'+fmt(sample.mu_hat-1, 6)]);
    addRow('#score-table tbody', [Number(size).toLocaleString('en-US'),
      fmt(sample.before.maximum_absolute_nominal_alpha_score, 6), fmt(corrected.largest_nominal_systematic_score, 6)]);
    addRow('#comparison-table tbody', [Number(size).toLocaleString('en-US'),
      fmt(sample.comparison.max_abs_delta_q, 4), fmt(sample.comparison.rms_delta_q, 4),
      fmt(sample.comparison.largest_residual.mu, 2), fmt(sample.rows.find(r => r.mu === .01).q, 4)]);

  }
  const small = data.samples['4096'], large = data.samples['8192'];
  const elapsed = seconds => `${Math.floor(seconds/60)}m ${String(Math.round(seconds%60)).padStart(2,'0')}s`;
  const phaseNames = {build_global_validate: 'Build + global fits + checks', scan: '4 parallel profile chunks', finalize: 'Merge + profile validation'};
  const phaseOrder = Object.keys(phaseNames);
  const resourceRows = [...data.resources.phases].sort((a,b) => phaseOrder.indexOf(a.phase)-phaseOrder.indexOf(b.phase) || a.reference_events-b.reference_events);
  for (const r of resourceRows) addRow('#resource-table tbody', [Number(r.reference_events).toLocaleString('en-US'), phaseNames[r.phase], elapsed(r.wall_seconds), fmt(r.total_cpu_seconds/60, 2), fmt(r.largest_job_peak_rss_gib, 2)+' GiB']);
  document.querySelector('#sample-difference').textContent = `Between the two corrected samples: maximum |Δq| = ${fmt(data.sample_difference.max_abs_delta_q, 4)} at μ = ${fmt(data.sample_difference.largest_residual.mu, 2)}; RMS Δq = ${fmt(data.sample_difference.rms_delta_q, 4)} across all 31 coordinates.`;
  for (const r of large.before.control_regions) addRow('#cr-table tbody', [r.channel, fmt(r.combined_data_count, 6), fmt(r.model_expected_count, 6), r.data_minus_model.toExponential(2)]);
  for (const [label, key] of [['Dataset rows (SR + 3 CR)', 'dataset_rows'], ['Total weighted count', 'sum_weights'], ['Negative event weights retained', 'negative_weights']]) {
    addRow('#construction-table tbody', [label, ...[small, large].map(s => key === 'sum_weights' ? fmt(s.construction[key], 6) : s.construction[key])]);
  }
  addRow('#construction-table tbody', ['Independent ΔNLL agreement', ...[small, large].map(s => s.construction.max_delta_nll_difference.toExponential(2))]);
  const alphaSmall = small.parameters.filter(p => p.name.startsWith('alpha_')).sort((a,b) => Math.abs(b.value) - Math.abs(a.value));
  const alphaLarge = Object.fromEntries(large.parameters.map(p => [p.name, p.value]));
  for (const p of alphaSmall.slice(0, 6)) addRow('#nuisance-table tbody', [p.name.replace('alpha_ATLAS_', ''), fmt(p.value, 4), fmt(alphaLarge[p.name], 4)]);
  document.querySelector('#fit-state-note').textContent = 'Both starts converged with status 0 and EDM below 6 × 10⁻⁸. Saved covariance estimates come from MIGRAD; no explicit full Hesse calculation was run.';
  document.querySelector('#snapshot-note').textContent = data.status === 'complete'
    ? 'Both scans and their numerical validation are complete in this snapshot.'
    : 'Asimov scans and final checks were still in progress at this snapshot. Counts and tables are not live; update the export to incorporate later results.';
  const commit = document.createElement('a');
  commit.href = `https://github.com/rafaellopesdesa/offshell_reweighting/commit/${data.analysis_commit}`;
  commit.textContent = `Analysis commit ${data.analysis_commit.slice(0, 7)}`;
  document.querySelector('#source-commit').appendChild(commit);
  const dialog = document.querySelector('#contents');
  const sections = [...document.querySelectorAll('.slides > section')];
  sections.forEach((section, index) => {
    const li = document.createElement('li'); const a = document.createElement('a');
    a.href = reading ? `#${section.id}` : `#/${section.id}`; a.textContent = section.querySelector('h1,h2').textContent;
    a.addEventListener('click', () => { dialog.close(); if (reading) section.scrollIntoView(); else Reveal.slide(index); });
    li.appendChild(a); document.querySelector('#contents-list').appendChild(li);
  });
  document.querySelector('#contents-button').addEventListener('click', () => dialog.showModal());
  document.querySelector('#close-contents').addEventListener('click', () => dialog.close());
  if (!reading) Reveal.initialize({ width: 1200, height: 720, margin: .055, minScale: .15, maxScale: 1.7,
    center: false, hash: true, slideNumber: 'c/t', controls: true, progress: true,
    transition: 'fade', transitionSpeed: 'fast', keyboardCondition: () => !dialog.open,
    postMessage: false, pdfSeparateFragments: false, view: 'slide' });
})();
