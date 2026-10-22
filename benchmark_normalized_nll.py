#!/usr/bin/env python3
"""Time fixed likelihood updates on one allocated CPU, excluding initialization."""
import argparse
import json
import math
import os
from pathlib import Path
import platform
import resource
import statistics
import time

from root_workspace_fit import load_model, write_json


def main(args):
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Run the benchmark in a Slurm allocation")
    cpu=min(os.sched_getaffinity(0))
    os.sched_setaffinity(0,{cpu})
    source,w,mc,data,parameters,floating,nll=load_model(args)
    values={r["name"]:r["value"] for r in json.loads(args.parameters.read_text())["parameters"]}
    values["mu"]=1.
    for name,value in values.items():w.var(name).setVal(value)
    alpha=sorted(n for n in values if n.startswith("alpha_"))
    def update(kind,i):
        delta=(1 if i%2 else -1)*1e-3
        if kind=="one_alpha":
            name=alpha[i%len(alpha)];w.var(name).setVal(values[name]+delta)
        elif kind=="all_alpha":
            for j,name in enumerate(alpha):w.var(name).setVal(values[name]+delta*math.sin(j+1))
        else:w.var("mu").setVal(1+delta)
    records={}
    for kind,count in [("one_alpha",128),("all_alpha",32),("mu_only",64)]:
        for name,value in values.items():w.var(name).setVal(value)
        for i in range(8):
            update(kind,i);assert math.isfinite(nll.getVal())
        wall,cpu_times=[],[]
        for i in range(count):
            update(kind,i+8)
            t=time.perf_counter();c=time.process_time();value=nll.getVal()
            cpu_times.append(time.process_time()-c);wall.append(time.perf_counter()-t)
            if not math.isfinite(value):raise RuntimeError("Invalid timed likelihood point")
        records[kind]={"evaluations":count,"median_wall_seconds":statistics.median(wall),
                       "mean_wall_seconds":statistics.mean(wall),"mean_cpu_seconds":statistics.mean(cpu_times),
                       "wall_seconds":wall}
    report={"label":args.label,"input":str(args.input.resolve()),"job_id":os.environ["SLURM_JOB_ID"],
            "hostname":platform.node(),"cpu_affinity":cpu,"peak_rss_gib":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/2**20,
            "records":records,"method":"Same CPU and trial parameter sequence; native CPU backend; first evaluations and initialization excluded. Data weights differ between the two workspaces, so this measures practical model evaluation cost, not a full profile or an isolated class microbenchmark."}
    write_json(args.output,report)
    print(json.dumps({"label":args.label,"records":{k:{n:v for n,v in r.items() if n!='wall_seconds'} for k,r in records.items()}}),flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input",type=Path,required=True)
    p.add_argument("--parameters",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--label",required=True)
    p.add_argument("--workspace",default="combined")
    p.add_argument("--model-config",default="ModelConfig_4l_offshell")
    p.add_argument("--data",default="obsData_4l_offshell")
    p.set_defaults(weighted_asimov=True)
    main(p.parse_args())
