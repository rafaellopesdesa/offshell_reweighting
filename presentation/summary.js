/* Frozen scientific results; no live cluster requests. */
(() => {
  'use strict';
  const data = window.SUMMARY_RESULTS;
  const family = window.FAMILY_RESULTS;
  if (!data?.complete || !family?.complete || family.models.length !== 6) {
    throw new Error('A complete validated summary is required.');
  }
  const reading = window.matchMedia('(max-width:700px)').matches;
  document.body.classList.toggle('reading', reading);
  const fmt = (n, digits = 2) => Number(n).toFixed(digits);
  const set = (selector, value, root = document) => {
    root.querySelectorAll(selector).forEach(el => { el.textContent = value; });
  };
  set('#snapshot-stamp', 'VALIDATED · ' + new Date(data.created_utc).toISOString().slice(0,16).replace('T',' ') + ' UTC');
  for (const m of family.models) {
    const section = document.querySelector(`[data-family-model="${m.label}"]`);
    set('[data-family-mu]', fmt(m.mu_hat, 7), section);
    set('[data-family-rms]', fmt(m.rms_delta_q, 5), section);
    set('[data-family-max]', fmt(m.max_abs_delta_q, 5), section);
    set('[data-family-max-mu]', fmt(m.mu_at_max_abs_delta_q), section);
  }
  for (const m of data.resources.models) {
    const section = document.querySelector(`[data-resource-model="${m.label}"]`);
    section.querySelectorAll('[data-resource]').forEach(el => { el.textContent = fmt(m[el.dataset.resource]); });
    const tr = document.createElement('tr');
    [m.title + ' · ' + m.reference_events.toLocaleString('en-US'),
      m.median_wall_minutes, m.minimizer_cpu_hours_including_overlaps,
      m.whole_job_cpu_hours, m.scan_wall_hours, m.maximum_slurm_peak_rss_gib].forEach(value => {
      const td = document.createElement('td'); td.textContent = typeof value === 'number' ? fmt(value) : value; tr.appendChild(td);
    });
    document.querySelector('#resource-summary-table tbody').appendChild(tr);
  }
  const medians = data.resources.models.map(m => m.median_wall_minutes);
  set('#median-cost-range', `${fmt(Math.min(...medians),1)}–${fmt(Math.max(...medians),1)}`);
  set('#maximum-process-memory', fmt(Math.max(...data.resources.models.map(m => m.maximum_process_peak_rss_gib))));
  for (const [size, sample] of Object.entries(data.diagnostics.samples)) {
    set(`[data-max-pull="${size}"]`, sample.maximum_absolute_pull_from_asimov.toExponential(2));
    set(`[data-hesse-sigma="${size}"]`, fmt(sample.hessian_sigma, 4));
    set(`[data-closure="${size}"]`, fmt(100*sample.closure_relative, 4) + '%');
    for (const kind of ['pulls','uncertainties']) {
      const section = document.querySelector(`#${kind}-${size}`);
      section.querySelectorAll('[data-diagnostic-view]').forEach(button => {
        button.addEventListener('click', () => {
          const view = button.dataset.diagnosticView;
          const suffix = kind === 'pulls' ? 'pulls_' + view : view;
          const base = `assets/summary_data_like_${size}_${suffix}`;
          section.querySelector('[data-diagnostic-figure]').src = base + '.svg';
          section.querySelector('[data-diagnostic-figure]').alt = `Data-like unpruned ${size}: ` + {
            nominal:'nuisance estimates and injected values in original nominal coordinates',
            asimov:'nuisance pulls relative to the Asimov auxiliary centers',
            auxiliary_impacts:'signed covariance responses to unit Gaussian auxiliary shifts',
            source_groups:'primary and grouped auxiliary uncertainty contributions'
          }[view];
          section.querySelector('[data-diagnostic-pdf]').href = base + '.pdf';
          set('[data-diagnostic-caption]', {
            nominal:'Original-nominal coordinates. Hollow diamonds mark the injected data-fit values.',
            asimov:'Asimov-centered pulls. Hollow diamonds are at zero; free qqZZ factors remain in raw units.',
            auxiliary_impacts:'Shared leading-source order at both sizes. Filled/outline bars show ±1σ auxiliary shifts.',
            source_groups:'All primary and auxiliary groups; source and Hesse totals are shown separately.'
          }[view], section);
          section.querySelectorAll('[data-diagnostic-view]').forEach(b => b.setAttribute('aria-pressed', String(b === button)));
        });
      });
    }
  }
  for (const figure of document.querySelectorAll('[data-paper-slot]')) {
    const config = (window.PAPER_FIGURES || {})[figure.dataset.paperSlot];
    if (!config) continue;
    const image = document.createElement('img');
    image.src = typeof config === 'string' ? config : config.src;
    image.alt = config.alt || 'Corresponding result from the paper';
    figure.querySelector('.paper-slot').replaceChildren(image);
    if (config.caption) figure.querySelector('figcaption').textContent = config.caption;
  }
  const dialog = document.querySelector('#contents');
  const sections = [...document.querySelectorAll('.slides > section')];
  sections.forEach((section,index) => {
    const li = document.createElement('li'); const link = document.createElement('a');
    link.href = reading ? '#' + section.id : '#/' + section.id;
    link.textContent = section.querySelector('h1,h2').textContent;
    link.addEventListener('click', () => { dialog.close(); if (reading) section.scrollIntoView(); else Reveal.slide(index); });
    li.appendChild(link); document.querySelector('#contents-list').appendChild(li);
  });
  document.querySelector('#contents-button').addEventListener('click', () => dialog.showModal());
  document.querySelector('#close-contents').addEventListener('click', () => dialog.close());
  if (reading) {
    document.querySelector('.masthead > a').href = '#overview';
    document.querySelectorAll('[data-jump]').forEach(a => { a.href = '#' + a.dataset.jump; });
    const scrollToHash = () => {
      const id = location.hash.replace(/^#\/?/,'');
      const section = document.getElementById(id);
      if (section) { history.replaceState(null,'','#'+id); section.scrollIntoView(); }
    };
    window.addEventListener('load',scrollToHash);
    window.addEventListener('hashchange',scrollToHash);
  } else {
    Reveal.initialize({width:1200,height:720,margin:.055,minScale:.15,maxScale:1.7,center:false,
      hash:true,slideNumber:'c/t',controls:true,progress:true,transition:'fade',transitionSpeed:'fast',
      keyboardCondition:() => !dialog.open,postMessage:false,pdfSeparateFragments:false,view:'slide'});
  }
})();
