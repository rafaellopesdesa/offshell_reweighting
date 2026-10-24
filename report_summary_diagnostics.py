#!/usr/bin/env python3
"""Updated data-like 4096/8192 nuisance and source-covariance figures."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path

import numpy as np
from report_uncertainty_decomposition import analyze, csv_write


def figures(samples, names, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'svg.fonttype':'none',
        'svg.hashsalt':'offshell-summary-diagnostics','axes.spines.top':False,'axes.spines.right':False})
    color='#c95639'
    def label(name):
        return name.removeprefix('alpha_ATLAS_').replace('SIMPLIFIED_UncorrUncertainty','uncorr').replace('CorrUncertainty','corr')
    nominal_limit=math.ceil(max(abs(r['value'])+r['posterior_sigma'] for s in samples.values()
                              for r in s['impacts'] if r['name'] in names)*2)/2
    centered_limit=max(1.3, math.ceil(10*max(abs(r['pull_from_anchor'])+r['posterior_sigma']
        for s in samples.values() for r in s['impacts'] if r['name'] in names))/10+.1)
    impact_limit=1.15*max(abs(r['unit_auxiliary_response']) for s in samples.values()
                         for r in s['impacts'] if r['name'] in names)
    total_limit=1.12*max(max(s['hessian_sigma'],s['propagated_sigma']) for s in samples.values())
    def save(fig,name):
        for ext in ('svg','pdf','png'):
            fig.savefig(output/f'{name}.{ext}',dpi=180,metadata={'Date':None} if ext=='svg' else None)
        plt.close(fig)
    for size,s in samples.items():
        lookup={r['name']:r for r in s['impacts']};rows=[lookup[n] for n in names]
        qq=[r for r in s['impacts'] if not r['constrained']]
        for center in (False,True):
            fig,(ax,qax)=plt.subplots(2,1,figsize=(7.3,6.5),gridspec_kw={'height_ratios':[5,1],'hspace':.37})
            x=[r['pull_from_anchor'] if center else r['value'] for r in rows]
            injection=[0. if center else r['anchor'] for r in rows]
            y=np.arange(len(rows))
            if center:ax.axvspan(-1,1,color='#edf2f5',zorder=0)
            ax.axvline(0,color='#8997a3',lw=.8)
            ax.errorbar(x,y,xerr=[r['posterior_sigma'] for r in rows],fmt='o',ms=3.5,color=color,capsize=2,lw=1.1)
            ax.scatter(injection,y,marker='D',facecolors='white',edgecolors='#253746',s=16,zorder=4)
            ax.set_yticks(y,[label(n) for n in names],fontsize=8.3);ax.invert_yaxis()
            ax.set_xlim((-centered_limit,centered_limit) if center else (-nominal_limit,nominal_limit))
            ax.set_xlabel('α̂ − αA  (Asimov constraint units)' if center else 'α̂  (original nominal units)')
            ax.grid(axis='x',alpha=.16)
            ax.set_title(f'Data-like, unpruned · {size:,} events',fontsize=12,pad=30)
            ax.legend(handles=[Line2D([],[],color=color,marker='o',ms=4,label='Fit ± Hesse σ'),
                               Line2D([],[],color='#253746',marker='D',mfc='white',ls='none',ms=4,label='Injected value')],
                      loc='lower center',bbox_to_anchor=(.5,1.005),frameon=False,ncols=2,fontsize=8.5)
            qy=np.arange(3)
            qax.errorbar([r['value'] for r in qq],qy,xerr=[r['posterior_sigma'] for r in qq],fmt='o',ms=4,color=color,capsize=2)
            qax.scatter([r['injected_value'] for r in qq],qy,marker='D',facecolors='white',edgecolors='#253746',s=16,zorder=4)
            qax.set_yticks(qy,['qqZZ β₀','qqZZ β₁','qqZZ β₂'],fontsize=9);qax.invert_yaxis()
            qax.set_xlim(.65,1.3);qax.set_xlabel('Free normalization factors; no Gaussian pull');qax.grid(axis='x',alpha=.16)
            fig.subplots_adjust(left=.37,right=.98,bottom=.12,top=.88)
            fig.text(.03,.025,'15 leading sources, shared order. Compact labels; all 130 full NP names and values in CSV.',fontsize=8,color='#536777')
            save(fig,f'data_like_{size}_pulls_'+('asimov' if center else 'nominal'))
        fig,ax=plt.subplots(figsize=(7.3,6.5));y=np.arange(len(rows));response=np.array([r['unit_auxiliary_response'] for r in rows])
        ax.barh(y,response,height=.55,color=color,alpha=.85)
        ax.barh(y,-response,height=.55,facecolor='none',edgecolor=color,lw=1.1)
        ax.axvline(0,color='#8997a3',lw=.8);ax.set_xlim(-impact_limit,impact_limit)
        ax.set_yticks(y,[label(n) for n in names],fontsize=8.3);ax.invert_yaxis();ax.grid(axis='x',alpha=.16)
        ax.set_xlabel('Local Δμ̂ from an auxiliary-observation shift')
        ax.set_title(f'Data-like, unpruned · {size:,} events',fontsize=12,pad=35)
        ax.legend(handles=[Patch(color=color,label='+1σ auxiliary'),Patch(facecolor='none',edgecolor=color,label='−1σ auxiliary')],
                  loc='lower center',bbox_to_anchor=(.5,1.005),frameon=False,ncols=2,fontsize=8.5)
        fig.subplots_adjust(left=.37,right=.98,bottom=.15,top=.88)
        fig.text(.035,.045,f"Primary σμ = {s['primary_sigma']:.4f}   Auxiliary σμ = {s['auxiliary_sigma']:.4f}   Hesse σμ = {s['hessian_sigma']:.4f}",fontsize=8.5)
        fig.text(.035,.02,'Linear covariance response V(μ,α); not finite nuisance refits. Source variances add in quadrature.',fontsize=8,color='#536777')
        save(fig,f'data_like_{size}_auxiliary_impacts')
        groups=[r for r in s['groups'] if r['kind'] in ('auxiliary','primary')]
        groups=sorted(groups,key=lambda r:r['variance'],reverse=True)
        totals=[r for r in s['groups'] if r['source'] in ('Propagated source total','Observed-Hessian total')]
        rows_group=groups+totals;y=np.arange(len(rows_group));fig,ax=plt.subplots(figsize=(7.3,6.5))
        colors=['#126c96' if r['kind']=='primary' else color if r['kind']=='auxiliary' else '#253746' for r in rows_group]
        ax.barh(y,[r['sigma'] for r in rows_group],height=.6,color=colors,alpha=.85)
        ax.axhline(len(groups)-.5,color='#9aa6af',lw=.8)
        ax.set_yticks(y,[r['source'] for r in rows_group],fontsize=9);ax.invert_yaxis();ax.set_xlim(0,total_limit);ax.grid(axis='x',alpha=.16)
        ax.set_xlabel('Contribution to σμ  (combine sources in quadrature)')
        ax.set_title(f'Data-like, unpruned · {size:,} events',fontsize=12,pad=22)
        ax.legend(handles=[Patch(color='#126c96',label='Primary data'),Patch(color=color,label='Gaussian auxiliaries')],
                  loc='lower center',bbox_to_anchor=(.5,1.0),frameon=False,ncols=2,fontsize=8.5)
        fig.subplots_adjust(left=.35,right=.98,bottom=.15,top=.86)
        fig.text(.035,.05,f"Source/Hesse variance difference: {100*s['closure_relative']:+.4f}%. No rescaling to enforce closure.",fontsize=8.5)
        note='Signed quadrature estimate: negative integration terms retained.' if s['negative_rows'] else 'Positive quadrature masses supply the generating Poisson measure.'
        fig.text(.035,.025,note+' Local errors, not profile intervals.',fontsize=8,color='#536777')
        save(fig,f'data_like_{size}_source_groups')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--results',type=Path,required=True)
    p.add_argument('--data-fit',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();assert os.environ.get('SLURM_JOB_ID'),'Use a Slurm allocation'
    args.output.mkdir(parents=True,exist_ok=True)
    observed=json.loads(args.data_fit.read_text());injected={r['name']:r['value'] for r in observed['parameters']}
    samples={}
    for size in (4096,8192):
        directory=args.results/f'data_like_{size}';s=analyze(directory)
        response=json.loads((directory/'observation_response.json').read_text())
        s['negative_rows']=response.get('quadrature',{}).get('negative_rows',0)
        s['quadrature']=response.get('quadrature',{'reference_events':size,'negative_rows':0})
        s['primary_variance_convention']=response['primary_variance_convention']
        for r in s['impacts']:
            r['injected_value']=injected[r['name']]
            if r['constrained']:assert abs(r['anchor']-r['injected_value'])<1e-12
        s['maximum_absolute_pull_from_asimov']=max(abs(r['pull_from_anchor']) for r in s['impacts'] if r['constrained'])
        s['reused_covariance']=(directory/'reuse.json').exists()
        samples[size]=s
    maps={size:{r['name']:r for r in s['impacts']} for size,s in samples.items()}
    constrained=[r['name'] for r in samples[4096]['impacts'] if r['constrained']]
    common=sorted(constrained,key=lambda n:max(abs(maps[k][n]['unit_auxiliary_response']) for k in samples),reverse=True)
    free=sorted(n for n,r in maps[4096].items() if not r['constrained'])
    ranking=[]
    for size,s in samples.items():
        ranks={name:i+1 for i,name in enumerate(sorted(constrained,key=lambda n:abs(maps[size][n]['unit_auxiliary_response']),reverse=True))}
        rows=[]
        for i,name in enumerate(common+free):
            r=maps[size][name]
            row={'reference_events':size,'shared_rank':i+1 if r['constrained'] else None,
                 'model_rank':ranks.get(name),'name':name,**{k:r[k] for k in ['constrained','group','value','injected_value','pull_from_anchor',
                 'posterior_sigma','unit_auxiliary_response','auxiliary_variance','conditional_postfit_np_impact','correlation_with_mu']}}
            rows.append(row);ranking.append(row)
        csv_write(args.output/f'data_like_{size}_all_nps.csv',rows)
        csv_write(args.output/f'data_like_{size}_uncertainty_sources.csv',s['groups'])
        cov=json.loads((args.results/f'data_like_{size}'/'covariance.json').read_text())
        csv_write(args.output/f'data_like_{size}_covariance.csv',[{'parameter':n,**dict(zip(cov['parameter_order'],v))} for n,v in zip(cov['parameter_order'],cov['covariance'])])
    csv_write(args.output/'all_np_ranking_4096_8192.csv',ranking)
    figures(samples,common[:15],args.output)
    summary={'complete':True,'created_utc':datetime.now(timezone.utc).isoformat(),'samples':samples,'displayed_nps':common[:15],
        'reference':'https://arxiv.org/abs/2608.15312','source_definition':'V[mu,alpha]^2 for each independent unit Gaussian auxiliary; primary score integral propagated with V[:,mu]',
        'pull_definition':'alpha_hat minus the actual Asimov auxiliary center; original-nominal-coordinate view separately shows injected data-fit shifts',
        'closure_policy':'Report source and observed-Hesse totals and directional checks separately; never rescale to enforce closure',
        'signed_quadrature_policy':'8192 primary source terms approximate the population score integral; negative integration weights are preserved and are not interpreted as negative Poisson variances',
        'observed_fit_sha256':hashlib.sha256(args.data_fit.read_bytes()).hexdigest()}
    (args.output/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    print(json.dumps({size:{k:s[k] for k in ['hessian_sigma','primary_sigma','auxiliary_sigma','closure_relative','negative_rows']} for size,s in samples.items()},indent=2))


if __name__=='__main__':main()
