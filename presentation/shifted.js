/* Validated data-like study; the older slides retain their historical snapshot. */
(() => {
  'use strict';
  const d = window.SHIFTED_RESULTS;
  if (!d || !d.complete) throw new Error('A complete data-like snapshot is required.');
  const fmt = (v, n=4) => Number(v).toFixed(n);
  const sci = v => Number(v).toExponential(2);
  const row = (id, values) => {
    const tr=document.createElement('tr');
    for (const value of values) { const td=document.createElement('td'); td.textContent=value; tr.appendChild(td); }
    document.querySelector(id+' tbody').appendChild(tr);
  };
  const names=['mu','mu_qqZZ','mu_qqZZ_1','mu_qqZZ_2'];
  for (const name of names) row('#shifted-injection-table',[name,fmt(d.generation.generation_parameters[name],7)]);
  row('#shifted-count-table',['SR',fmt(d.generation.sr_asimov_yield,3)]);
  d.generation.cr_asimov_counts.forEach((v,i)=>row('#shifted-count-table',['CR'+i,fmt(v,3)]));
  const fit=d.samples.data_like;
  for (const a of fit.attempts) row('#shifted-fit-table',[fmt(a.start_mu,1),fmt(a.mu_hat,8),sci(a.edm),'Quality '+a.covariance_quality+' · '+a.role]);
  document.querySelector('#shifted-score').textContent=sci(d.construction_validation.maximum_absolute_score);
  document.querySelector('#shifted-closure').textContent=sci(d.construction_validation.maximum_process_integral_error);
  document.querySelector('#shifted-recovery').textContent=`Selected μ̂ = ${fmt(fit.mu_hat,8)}. Largest parameter displacement from injection: ${sci(d.recovery.maximum_parameter_difference)}; NLL − NLL at generation: ${sci(d.recovery.nll_minus_generation)}.`;
  for (const [key,label] of [['true','True'],['data_like','Data-like']]) {
    const s=d.samples[key];
    row('#shifted-comparison-table',[label,fmt(s.mu_hat,6),fmt(s.comparison.full.rms_delta_q,4),fmt(s.q_zero,4),fmt(s.q_two,4)]);
  }
  document.querySelector('#shifted-likelihood-check').textContent=`All 42 points in each scan converge. Maximum independent ΔNLL disagreement: ${sci(Math.max(...Object.values(d.samples).map(s=>s.maximum_native_explicit_disagreement)))}. Data-like SR closure over the scan: ${sci(fit.maximum_sr_integral_mismatch)}.`;
  const phaseNames={upper_extension:'Extend true to μ = 2',global_workflow:'Build + global fits + checks',scan:'5 parallel profile chunks',finalize:'Merge + profile validation'};
  const duration=s=>`${Math.floor(s/60)}m ${String(Math.round(s%60)).padStart(2,'0')}s`;
  for (const r of d.resources.phases) row('#shifted-resource-table',[r.sample==='true'?'True':'Data-like',phaseNames[r.phase],duration(r.wall_seconds),fmt(r.total_cpu_seconds/60,2),fmt(r.largest_job_peak_rss_gib,2)+' GiB']);
  for (const [key,label] of [['one_alpha','One α changes'],['all_alpha','All 127 α change'],['mu_only','Only μ changes']]) {
    const r=d.timing.updates[key];
    row('#normalizer-evaluation-table',[label,fmt(1000*r.true,3),fmt(1000*r.data_like,3),fmt(r.data_like_over_true,3)+'×']);
  }
  for (const [key,label] of [['true','True'],['data_like','Data-like']]) {
    const r=d.timing.profiles[key];
    row('#normalizer-profile-table',[label,duration(r.seconds),r.likelihood_evaluations,sci(r.edm)]);
  }
  document.querySelector('#normalizer-performance-note').textContent=`On ${d.timing.hostname}, CPU ${d.timing.cpu_affinity}: complete μ = 0.5 profile time ratio ${fmt(d.timing.profile_time_ratio,3)}×; likelihood-call ratio ${fmt(d.timing.profile_call_ratio,3)}×. Evaluation-benchmark peak RSS: ${fmt(d.timing.peak_rss_gib.true,2)} → ${fmt(d.timing.peak_rss_gib.data_like,2)} GiB.`;
  document.querySelector('#snapshot-stamp').textContent='COMPLETE · '+d.created_utc.replace('T',' ').replace('+00:00',' UTC');
})();
