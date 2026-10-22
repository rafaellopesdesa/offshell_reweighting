#!/usr/bin/env python3
"""Check saved Hesse curvature with independent native-NLL directional differences."""
import argparse
import json
from pathlib import Path

import numpy as np
import ROOT

from root_workspace_fit import load_model,restore,sha256,write_json


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input",type=Path,required=True)
    p.add_argument("--bestfit",type=Path,required=True)
    p.add_argument("--covariance",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--workspace",default="combined")
    p.add_argument("--model-config",default="ModelConfig_4l_offshell")
    p.add_argument("--data",default="obsData_4l_offshell")
    p.add_argument("--weighted-asimov",action="store_true")
    args=p.parse_args()
    source,w,mc,data,parameters,floating,nll=load_model(args)
    best=json.loads(args.bestfit.read_text());cov=json.loads(args.covariance.read_text())
    assert cov["passed"] and best["input_sha256"]==cov["input_sha256"]==sha256(args.input)
    restore(w,best["parameters"])
    names=cov["parameter_order"];point={r["name"]:r["value"] for r in cov["parameters"]}
    values=np.array([point[n] for n in names]);V=np.array(cov["covariance"])
    def evaluate(p):
        for name,value in zip(names,p):w.var(name).setVal(float(value))
        return nll.getVal()
    baseline=evaluate(values)
    im=names.index("mu");eigen,rot=np.linalg.eigh(V);root=rot@np.diag(np.sqrt(eigen))
    rng=np.random.default_rng(260815312)
    directions={"profile_mu_tangent":V[:,im]/np.sqrt(V[im,im]),
                "largest_covariance_eigenvector":rot[:,-1]*np.sqrt(eigen[-1])}
    for i in range(2):
        unit=rng.normal(size=len(names));unit/=np.linalg.norm(unit)
        directions[f"whitened_direction_{i}"]=root@unit
    checks=[]
    for name,direction in directions.items():
        differences=[]
        for h in (.02,.01):
            plus=evaluate(values+h*direction);minus=evaluate(values-h*direction)
            differences.append((plus+minus-2*baseline)/h**2)
        extrap=(4*differences[1]-differences[0])/3
        checks.append({"direction":name,"h_002":differences[0],"h_001":differences[1],
                       "richardson_curvature":extrap,"expected_curvature":1.,"disagreement":extrap-1})
    evaluate(values)
    result={"passed":all(abs(c["disagreement"])<.005 for c in checks),
            "covariance_sha256":sha256(args.covariance),"checks":checks,
            "maximum_absolute_disagreement":max(abs(c["disagreement"]) for c in checks)}
    write_json(args.output,result);source.Close()
    print(json.dumps(result,indent=2))
    if not result["passed"]:raise RuntimeError("Directional Hesse validation failed")


if __name__=="__main__":main()
