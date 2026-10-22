/* Additive physical sources, with conditional nuisance diagnostics kept separate. */
(() => {
  'use strict';
  const d=window.DECOMPOSITION_RESULTS;
  if (!d || !d.complete || !d.resources.complete) throw new Error('Complete covariance analysis required.');
  const fmt=(v,n=4)=>Number(v).toFixed(n);
  const row=(id,values)=>{
    const tr=document.createElement('tr');
    for(const value of values){const td=document.createElement('td');td.textContent=value;tr.appendChild(td);}
    document.querySelector(id+' tbody').appendChild(tr);
  };
  const keys=['true','data_like'];
  const values=name=>keys.map(k=>d.samples[k].groups.find(r=>r.source===name));
  const primary=['SR statistics','CR0 statistics','CR1 statistics','CR2 statistics','Primary-data subtotal','Auxiliary subtotal','Propagated source total','Observed-Hessian total'];
  for(const name of primary)row('#uncertainty-primary',[name,...values(name).map(r=>fmt(r.sigma))]);
  const groups=d.samples.data_like.groups.filter(r=>r.kind==='auxiliary').map(r=>r.source);
  for(const name of groups)row('#uncertainty-auxiliary',[name,...values(name).map(r=>fmt(r.sigma))]);
  row('#uncertainty-closure',['Source-sum / Hessian variance − 1',...keys.map(k=>fmt(100*d.samples[k].closure_relative,3)+'%')]);
  row('#uncertainty-closure',['Largest curvature-check discrepancy',...keys.map(k=>fmt(100*d.samples[k].curvature_validation.maximum_absolute_disagreement,3)+'%')]);
  row('#uncertainty-closure',['Source / native-curvature probe − 1',...keys.map(k=>fmt(100*d.samples[k].source_to_native_curvature_probe_relative_difference,4)+'%')]);
  row('#uncertainty-closure',['qqZZ block freeze/float Δvariance',...keys.map(k=>fmt(d.samples[k].conditional_diagnostics.qqzz_block_profiling_variance_penalty))]);
  const duration=s=>`${Math.floor(s/60)}m ${String(Math.round(s%60)).padStart(2,'0')}s`;
  for(const r of d.resources.totals)row('#uncertainty-resources',[r.model==='true'?'True':'Data-like',duration(r.wall_seconds),fmt(r.cpu_seconds/60,2),fmt(r.largest_job_peak_rss_gib,2)+' GiB']);
  document.querySelector('#uncertainty-check-note').textContent=`The source/native ratio uses curvature along V[:,μ]. No sources are rescaled. Four covariance directions pass two-scale checks; maximum primary-variance step-size change: ${Number(Math.max(...keys.map(k=>d.samples[k].maximum_step_relative_difference))).toExponential(2)}.`;
})();
