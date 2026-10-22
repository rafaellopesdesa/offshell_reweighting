#!/usr/bin/env python3
"""Summarize same-CPU evaluation and single-profile timings."""
import argparse
import hashlib
import json
from pathlib import Path
import statistics


def main(args):
    files=sorted(args.results.glob('[0-3]_*.json'))
    if len(files)!=4:raise RuntimeError('Require alternating true/data-like/data-like/true runs')
    runs=[json.loads(p.read_text()) for p in files]
    assert len({(r['hostname'],r['cpu_affinity'],r['job_id']) for r in runs})==1
    assert [r['label'] for r in runs]==['true','data_like','data_like','true']
    report={'complete':True,'hostname':runs[0]['hostname'],'cpu_affinity':runs[0]['cpu_affinity'],
            'slurm_job_id':runs[0]['job_id'],'method':runs[0]['method'],'updates':{},'profiles':{},
            'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files}}
    for kind in runs[0]['records']:
        row={label:statistics.mean(r['records'][kind]['mean_cpu_seconds'] for r in runs if r['label']==label)
             for label in ['true','data_like']}
        row['data_like_over_true']=row['data_like']/row['true']
        report['updates'][kind]=row
    report['peak_rss_gib']={label:max(r['peak_rss_gib'] for r in runs if r['label']==label) for label in ['true','data_like']}
    for label in ['true','data_like']:
        path=args.results/f'profile_{label}'
        check=json.loads((path/'validation.json').read_text())
        assert check['passed']
        rows=json.loads((path/'scan_points.json').read_text())
        assert len(rows)==1 and rows[0]['mu']==.5
        row=rows[0]
        assert row['migrad_status']==0 and row['edm']<1e-4 and row['likelihood_evaluations']>0
        report['profiles'][label]={k:row[k] for k in ['mu','q','seconds','likelihood_evaluations','edm','migrad_status']}
        report['profiles'][label]['seconds_per_evaluation']=row['seconds']/row['likelihood_evaluations']
    report['profile_time_ratio']=report['profiles']['data_like']['seconds']/report['profiles']['true']['seconds']
    report['profile_call_ratio']=report['profiles']['data_like']['likelihood_evaluations']/report['profiles']['true']['likelihood_evaluations']
    report['profile_method']='One mu=0.5 point from each global-fit nuisance vector, 130 free nuisances, same Minuit settings and allocated CPU. Different Asimov data/model; call counts and profile curvatures can differ. Timing includes minimization only; loading and independent checks are excluded.'
    (args.results/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--results',type=Path,required=True)
    main(p.parse_args())
