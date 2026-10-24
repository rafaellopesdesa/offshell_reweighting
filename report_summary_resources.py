#!/usr/bin/env python3
"""Measured per-point and whole-job resources for the six validated scans."""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess

import numpy as np
from collect_anchor_resources import seconds, memory_gib


def write_csv(path, rows):
    with path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()



def write_public_report(report, computed, output):
    """Publish numerical costs with anonymous workers; retain cluster IDs locally."""
    import copy
    output.mkdir(parents=True, exist_ok=True)
    worker_ids={job['job_id']:i+1 for i,job in enumerate(report['jobs'])}
    assignments={row['slurm_job_id']:(row['workspace'],row['chunk']) for row in computed}
    public=copy.deepcopy(report)
    public.pop('accounting_command')
    public['accounting_method']='Slurm sacct; detailed job IDs, hosts and execution timestamps retained only in the local audit'
    public['public_export_source_sha256']=digest(Path(__file__))
    public['jobs']=[{'worker_index':worker_ids[job['job_id']],
        'workspace':assignments[job['job_id']][0],'profile_chunk':assignments[job['job_id']][1],
        **{key:job[key] for key in ('wall_seconds','cpu_seconds','slurm_peak_rss_gib','allocated_cpus')}}
        for job in report['jobs']]
    for model in public['models']:
        for row in model['rows']:
            row['worker_index']=worker_ids[row.pop('slurm_job_id')]
            row['chunk_index']=row.pop('array_task_id')
        write_csv(output/f"{model['label']}_resources.csv",model['rows'])
    public_computed=[{**{key:value for key,value in row.items() if key not in ('slurm_job_id','array_task_id')},
        'worker_index':worker_ids[row['slurm_job_id']]} for row in computed]
    write_csv(output/'all_computed_profiles.csv',public_computed)
    write_csv(output/'slurm_jobs.csv',public['jobs'])
    write_csv(output/'resource_summary.csv',[{k:v for k,v in m.items() if k not in ('rows','color')} for m in public['models']])
    (output/'summary.json').write_text(json.dumps(public,indent=2,allow_nan=False)+'\n')


def plots(models, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'svg.fonttype':'none',
        'svg.hashsalt':'offshell-summary-resources','axes.spines.top':False,'axes.spines.right':False})
    ymax = math.ceil(max(m['maximum_wall_minutes'] for m in models)/5)*5
    memory_max = math.ceil(max(m['maximum_process_peak_rss_gib'] for m in models)*2)/2+.25
    for model in models:
        rows=model['rows']; x=[r['mu'] for r in rows]
        fig=plt.figure(figsize=(11.5,6.3)); grid=fig.add_gridspec(2,2,height_ratios=[1.4,1],hspace=.5,wspace=.28)
        time=fig.add_subplot(grid[0,:]); cpu=fig.add_subplot(grid[1,0]); mem=fig.add_subplot(grid[1,1])
        time.plot(x,[r['wall_seconds']/60 for r in rows],'o-',ms=3,color=model['color'],label='Wall time')
        time.plot(x,[r['cpu_seconds']/60 for r in rows],'--',lw=1.6,color='#253746',label='CPU time')
        time.set_ylabel('Minimization minutes');time.set_ylim(0,ymax);time.legend(frameon=False,ncols=2)
        time.set_title(f"{model['title']} · {model['reference_events']:,} events · one CPU per profile",fontsize=13)
        cpu.plot(x,[r['cpu_utilization_percent'] for r in rows],'.-',color=model['color'],ms=4)
        cpu.set_ylabel('CPU / wall time (%)');cpu.set_ylim(0,110)
        mem.plot(x,[r['process_peak_rss_gib'] for r in rows],'.-',color=model['color'],ms=4)
        mem.set_ylabel('Process peak RSS (GiB)');mem.set_ylim(0,memory_max)
        for ax in (time,cpu,mem):
            ax.grid(alpha=.18);ax.set_xlim(0,3);ax.set_xticks(np.arange(0,3.01,.5));ax.set_xlabel('Signal strength μ')
        fig.subplots_adjust(left=.09,right=.985,bottom=.15,top=.92)
        fig.text(.09,.03,'Warm starts within chunks; initialization excluded from point times. RSS is a process high-water mark.',fontsize=9,color='#536777')
        for ext in ('svg','pdf','png'):
            fig.savefig(output/f"{model['label']}_resources.{ext}",dpi=180,metadata={'Date':None} if ext=='svg' else None)
        plt.close(fig)
    fig,(wall,mem)=plt.subplots(1,2,figsize=(11.5,5.2))
    for model in models:
        rows=model['rows'];x=[r['mu'] for r in rows]
        style='-' if model['reference_events']==4096 else '--'
        title=model['title'].replace(' Asimov','')+f" · {model['reference_events']}"
        wall.plot(x,[r['wall_seconds']/3600 for r in rows],style,color=model['color'],lw=1.5,label=title)
        mem.plot(x,[r['process_peak_rss_gib'] for r in rows],style,color=model['color'],lw=1.5)
    wall.axhline(24,color='#253746',ls=':',lw=1.8);wall.text(.06,27,'Paper: ≈24 h per point',fontsize=10)
    mem.axhline(500e9/2**30,color='#253746',ls=':',lw=1.8);mem.text(.06,540,'Paper: ≈500 GB (466 GiB)',fontsize=10)
    wall.set_ylabel('Minimization wall time (hours)');mem.set_ylabel('Process peak RSS (GiB)')
    wall.set_yscale('log');mem.set_yscale('log');wall.set_ylim(.001,60);mem.set_ylim(.7,900)
    wall.legend(frameon=False,fontsize=8.5,ncols=2,loc='lower right')
    for ax in (wall,mem):ax.set_xlabel('Signal strength μ');ax.set_xlim(0,3);ax.grid(alpha=.18)
    fig.suptitle('Per-point cost: current measurements and reported paper production',fontsize=13)
    fig.subplots_adjust(left=.08,right=.98,bottom=.2,top=.88,wspace=.28)
    fig.text(.08,.055,'Historical values supplied by the author. CPU count, hardware and RAM accounting are unspecified.\nCurrent profiles use one CPU; point times exclude loading. This is a contextual, not hardware-matched, comparison.',fontsize=9,color='#536777')
    for ext in ('svg','pdf','png'):fig.savefig(output/f'historical_resource_comparison.{ext}',dpi=180,metadata={'Date':None} if ext=='svg' else None)
    plt.close(fig)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    assert os.environ.get('SLURM_JOB_ID'),'Use a Slurm allocation'
    args.output.mkdir(parents=True,exist_ok=True)
    family=json.loads((args.run/'presentation_report/summary.json').read_text());assert family['complete']
    all_points=[];model_jobs={}
    for model in family['models']:
        rows=[]
        for path in sorted((args.run/model['label']/'scan').glob('chunk_*/scan_points.json')):
            for r in json.loads(path.read_text()):
                assert r['migrad_status']==0 and r['edm']<=1e-4
                rows.append({'workspace':model['label'],'chunk':path.parent.name,'mu':r['mu'],
                    'slurm_job_id':r['slurm_job_id'],'array_task_id':r['slurm_array_task_id'],
                    'wall_seconds':r['seconds'],'cpu_seconds':r['cpu_seconds'],
                    'process_peak_rss_gib':r['process_peak_rss_gib'],'likelihood_evaluations':r['likelihood_evaluations']})
        assert len(rows)==66
        all_points.extend(rows);model_jobs[model['label']]=sorted({r['slurm_job_id'] for r in rows})
        assert len(model_jobs[model['label']])==6
    fields=['JobID','JobIDRaw','State','ElapsedRaw','TotalCPU','MaxRSS','AllocCPUS','NodeList','Start','End','ExitCode']
    ids=sorted({r['slurm_job_id'] for r in all_points})
    command=['sacct','-n','-P','--units=K','-j',','.join(ids),'--format='+','.join(fields)]
    raw=subprocess.check_output(command,text=True,env={**os.environ,'TZ':'UTC'})
    (args.output/'slurm_accounting.txt').write_text('|'.join(fields)+'\n'+raw)
    records=[dict(zip(fields,line.split('|'))) for line in raw.splitlines() if line.strip()]
    lookup={r['JobIDRaw']:r for r in records};jobs=[]
    for id_ in ids:
        row=lookup[id_];assert row['State']=='COMPLETED' and row['ExitCode']=='0:0' and int(row['AllocCPUS'])==1
        peaks=[memory_gib(r['MaxRSS']) for r in records if r['JobIDRaw'].startswith(id_+'.') and r['MaxRSS']]
        assert peaks and seconds(row['TotalCPU'])>0
        jobs.append({'job_id':id_,'array_job_id':row['JobID'],'wall_seconds':int(row['ElapsedRaw']),
            'cpu_seconds':seconds(row['TotalCPU']),'slurm_peak_rss_gib':max(peaks),'allocated_cpus':1,
            'node':row['NodeList'],'start':row['Start'],'end':row['End']})
    models=[]
    for spec in family['models']:
        csvpath=args.run/'csv'/f"{spec['label']}.csv";assert digest(csvpath)==spec['scan_csv_sha256']
        selected=json.loads((args.run/spec['label']/'final/scan_points.json').read_text())
        lookup_points={r['mu']:r for r in selected};points=[]
        for rawrow in csv.DictReader(csvpath.open()):
            mu=float(rawrow['mu']);source=lookup_points[mu]
            row={'mu':mu,'wall_seconds':float(rawrow['wall_seconds']),'cpu_seconds':float(rawrow['cpu_seconds']),
                'process_peak_rss_gib':float(rawrow['process_peak_rss_gib']),
                'cpu_utilization_percent':100*float(rawrow['cpu_seconds'])/float(rawrow['wall_seconds']),
                'slurm_job_id':source['slurm_job_id'],'array_task_id':source['slurm_array_task_id']}
            assert abs(row['wall_seconds']-source['seconds'])<1e-10 and abs(row['cpu_seconds']-source['cpu_seconds'])<1e-10
            assert all(math.isfinite(row[k]) and row[k]>0 for k in ['wall_seconds','cpu_seconds','process_peak_rss_gib'])
            points.append(row)
        assert [r['mu'] for r in points]==[i/20 for i in range(61)]
        job_rows=[j for j in jobs if j['job_id'] in model_jobs[spec['label']]]
        computed=[r for r in all_points if r['workspace']==spec['label']]
        wall=[r['wall_seconds']/60 for r in points]
        model={k:spec[k] for k in ['label','title','reference_events','constrained_nps','color']}
        model.update(points=61,computed_points=66,median_wall_minutes=float(np.median(wall)),
            minimum_wall_minutes=min(wall),maximum_wall_minutes=max(wall),
            median_cpu_minutes=float(np.median([r['cpu_seconds']/60 for r in points])),
            maximum_process_peak_rss_gib=max(r['process_peak_rss_gib'] for r in points),
            maximum_slurm_peak_rss_gib=max(r['slurm_peak_rss_gib'] for r in job_rows),
            minimizer_cpu_hours_including_overlaps=sum(r['cpu_seconds'] for r in computed)/3600,
            whole_job_cpu_hours=sum(r['cpu_seconds'] for r in job_rows)/3600,
            scan_wall_hours=(max(datetime.fromisoformat(j['end']) for j in job_rows)-min(datetime.fromisoformat(j['start']) for j in job_rows)).total_seconds()/3600,
            scan_csv_sha256=spec['scan_csv_sha256'],rows=points)
        write_csv(args.output/f"{spec['label']}_resources.csv",points);models.append(model)
    report={'complete':True,'created_utc':datetime.now(timezone.utc).isoformat(),'models':models,'jobs':jobs,
        'historical_benchmark':{'wall_hours_per_point':24,'memory_gb':500,'memory_gib':500e9/2**30,
            'source':'Author report, confirmed 2026-10-24: one profiled mu point',
            'cpu_count':None,'hardware':None,'ram_accounting':'Unspecified: historical requirement, not certified peak RSS'},
        'definitions':{'point_time':'Minimization wall/CPU time, including strategy retry; excludes loading and initial NLL setup',
            'cpu_utilization':'100 * process CPU seconds / wall seconds; one allocated CPU per profile',
            'point_memory':'Process peak RSS up through the selected point within its warm-start chunk; not isolated per-point allocation',
            'aggregate_work':'All 66 profiles in six overlapping chunks per model; 61 unique coordinates shown',
            'job_memory':'Largest Slurm step MaxRSS per model, not sum of concurrently running processes',
            'scan_wall':'First chunk start to last chunk finish, including scheduling gaps; up to four chunks per model concurrently',
            'scope':'Profile scan only; global fits, final likelihood validation and separate Hesse calculations excluded',
            'comparison':'Historical and current hardware/start policy/accounting differ; no controlled hardware speedup claim'},
        'accounting_command':command,'source_sha256':digest(Path(__file__))}
    (args.output/'summary.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    write_csv(args.output/'resource_summary.csv',[{k:v for k,v in m.items() if k not in ['rows','color']} for m in models])
    write_csv(args.output/'all_computed_profiles.csv',all_points);write_csv(args.output/'slurm_jobs.csv',jobs)
    write_public_report(report,all_points,args.output/'public')
    plots(models,args.output)
    print(json.dumps([{k:m[k] for k in ['label','median_wall_minutes','maximum_wall_minutes','maximum_process_peak_rss_gib','maximum_slurm_peak_rss_gib']} for m in models],indent=2))


if __name__=='__main__':main()
