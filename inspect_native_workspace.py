#!/usr/bin/env python3
"""Read a native NRE workspace and export aggregate graph structure, never events."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import ROOT


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();ROOT.gROOT.SetBatch(True)
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
    source=ROOT.TFile.Open(str(args.input),'READ');w=source.Get('combined');assert w
    mc=w.obj('ModelConfig_4l_offshell');data=w.data('obsData_4l_offshell');assert mc and data
    nodes={}
    wanted=['simPdf','pdf_SR_4l_offshell','pdfSRUnc','pdf_process_S','prod_process_S',
        's_factor_S','s_factor_S_nr','s_factor_S_dr','piecewise_interp_S','mod_S',
        'nfFactor_S','nfFactor_SBI','nfFactor_B','nfFactor_EWB','nfFactor_EWSBI','nfFactor_EWSBI10',
        'nfFactor_qqZZ_0','nfFactor_qqZZ_1','nfFactor_qqZZ_2','nfFactor_ttV',
        'nominal_yield_SR_S','constraintPdf','gauss_ATLAS_EG_SCALE_ALL','pdf_CR0','constrainedCR0','nu_CR0_qqZZ']
    for name in wanted:
        node=w.arg(name);assert node,name
        item={'class':node.ClassName(),'servers':[v.GetName() for v in node.servers()]}
        item['server_count']=len(item['servers'])
        if len(item['servers'])>24:item['servers']=item['servers'][:6];item['servers_truncated']=True
        if node.InheritsFrom('RooFormulaVar'):item['formula']=node.GetTitle()
        if node.InheritsFrom('RooWrapperPdf'):item['self_normalized']=bool(node.selfNormalized())
        if name in ('mod_S','piecewise_interp_S'):item['interpolation_codes']=sorted(set(node.interpolationCodes()))
        nodes[name]=item
    channels={}
    for i in range(data.numEntries()):
        label=data.get(i).getCatLabel('cat');row=channels.setdefault(label,{'rows':0,'count':0.})
        row['rows']+=1;row['count']+=data.weight()
    sets={}
    for key,getter in [('observables',mc.GetObservables),('parameters_of_interest',mc.GetParametersOfInterest),
        ('nuisances',mc.GetNuisanceParameters),('global_observables',mc.GetGlobalObservables),
        ('conditional_observables',mc.GetConditionalObservables),('external_constraints',mc.GetExternalConstraints)]:
        values=getter();names=[v.GetName() for v in values] if values else []
        sets[key]={'count':len(names)}
        if len(names)<=16:sets[key]['names']=names
    globals_=[v.getVal() for v in mc.GetGlobalObservables()]
    assert all(v.isConstant() for v in mc.GetGlobalObservables())
    assert all(value==0 for value in globals_), 'The explanatory data example assumes original nominal centers'
    formulas={v.GetName():{'formula':v.GetTitle(),'servers':[x.GetName() for x in v.servers()]}
              for v in w.allFunctions() if v.InheritsFrom('RooFormulaVar')}
    with args.input.open('rb') as f:sha=hashlib.file_digest(f,'sha256').hexdigest()
    report={'source_basename':args.input.name,'source_sha256':sha,'root_version':ROOT.gROOT.GetVersion(),
        'workspace':w.GetName(),'dataset':data.GetName(),'data_rows':data.numEntries(),'sum_weights':data.sumEntries(),
        'channels':channels,'model_config_sets':sets,'pdf_classes':dict(Counter(v.ClassName() for v in w.allPdfs())),
        'function_classes':dict(Counter(v.ClassName() for v in w.allFunctions())),
        'selected_nodes':nodes,'formula_nodes':formulas,'global_observations_all_zero':True,'dataset_columns':data.get().getSize(),'scope':'Aggregate graph inspection of the supplied file; not the original production builder and no event values exported'}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('selected_nodes','formula_nodes')},indent=2))
    source.Close()


if __name__=='__main__':main()
