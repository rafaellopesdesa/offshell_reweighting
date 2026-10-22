#!/usr/bin/env python3
"""Inspect serialized normalizers and their connections to every SR component."""
import argparse
import json
import os
from pathlib import Path

import ROOT

from root_asimov_inputs import SR_PROCESSES
from root_workspace_fit import sha256, write_json
from same_sample_normalizer import load_normalizer


def main(args):
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Workspace audits require a Slurm allocation")
    ROOT.gROOT.SetBatch(True)
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
    load_normalizer()
    records = []
    for path in args.inputs:
        before = sha256(path)
        source = ROOT.TFile.Open(str(path), "READ")
        if not source or source.IsZombie():
            raise RuntimeError(f"Cannot read {path}")
        w = source.Get("combined")
        metadata = json.loads(str(source.Get("asimov_metadata").GetString()))
        mc = w.obj("ModelConfig_4l_offshell")
        data = w.data("obsData_4l_offshell")
        instances = sorted(str(f.GetName()) for f in w.allFunctions()
                           if f.ClassName() == "RooSameSampleNormalizer")
        connected = {p: bool(w.function(f"reference_Z_{p}") and
                            w.function(f"reference_Z_{p}").ClassName() == "RooSameSampleNormalizer" and
                            w.function(f"prod_process_{p}").dependsOn(w.function(f"reference_Z_{p}")))
                     for p in SR_PROCESSES}
        passed = len(instances) == len(SR_PROCESSES) and all(connected.values())
        record = dict(path=str(path.resolve()), sha256=before,
                      normalizer_instances=instances, connected_to_shapes=connected,
                      uses_full_normalization=passed,
                      metadata_implementation=metadata.get("normalization_implementation"),
                      reference_events=metadata["reference_events"],
                      rows=data.numEntries(), constrained_nps=mc.GetNuisanceParameters().getSize(),
                      floating_parameters=sum(not v.isConstant() for v in mc.GetPdf().getParameters(data)
                                              if v.InheritsFrom("RooRealVar")),
                      root_version=ROOT.gROOT.GetVersion())
        if args.nll:
            w.loadSnapshot("asimov_generation")
            nll = mc.GetPdf().createNLL(data, ROOT.RooFit.Extended(True),
                ROOT.RooFit.Constrain(mc.GetNuisanceParameters()),
                ROOT.RooFit.GlobalObservables(mc.GetGlobalObservables()),
                ROOT.RooFit.GlobalObservablesSource("model"),
                ROOT.RooFit.EvalBackend("cpu"), ROOT.RooFit.Offset("initial"))
            record["generation_nll"] = nll.getVal()
            record["generation_mu"] = w.var("mu").getVal()
            del nll
        source.Close()
        assert sha256(path) == before
        records.append(record)
        write_json(args.output, records)
        print(json.dumps(record), flush=True)
        if args.require_normalizers and not passed:
            raise RuntimeError(f"Missing or disconnected normalizers: {path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", type=Path, nargs="+")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--require-normalizers", action="store_true")
    parser.add_argument("--nll", action="store_true", help="Also evaluate the generating-point NLL")
    main(parser.parse_args())
