"""Small, resumable native RooFit scan for the three demonstration workspaces."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import resource
import time

import numpy as np
import ROOT


def sha256(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def save(path, value):
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def state(variables):
    return {v.GetName(): [v.getVal(), v.getError()] for v in variables}


def restore(variables, saved):
    for v in variables:
        v.setVal(saved[v.GetName()][0])
        v.setError(saved[v.GetName()][1])


def checked(row):
    return (row["status"] == 0 and math.isfinite(row["nll"])
            and math.isfinite(row["edm"]) and row["edm"] <= 1e-4)


def make_minimizer(nll, eps):
    minim = ROOT.RooMinimizer(nll)
    minim.setMinimizerType("Minuit2")
    minim.setEps(eps)
    minim.setErrorLevel(0.5)
    minim.setMaxFunctionCalls(100000)
    minim.setMaxIterations(10000)
    minim.setPrintLevel(-1)
    minim.setPrintEvalErrors(3)
    minim.setOffsetting(True)
    return minim


def minimize(nll, eps, minim=None):
    if minim is None:
        minim = make_minimizer(nll, eps)
    start, cpu = time.monotonic(), time.process_time()
    before = minim.evalCounter()
    for strategy in (0, 1):
        minim.setStrategy(strategy)
        status = minim.minimize("Minuit2", "Migrad")
        result = minim.save()
        row = dict(status=int(status), nll=float(nll.getVal()),
                   edm=float(result.edm()), covariance_quality=int(result.covQual()),
                   strategy=strategy)
        if checked(row):
            break
    row.update(seconds=time.monotonic() - start, cpu_seconds=time.process_time() - cpu,
               calls=int(minim.evalCounter() - before),
               process_peak_rss_gib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**2)
    if not checked(row):
        raise RuntimeError(f"Fit failed (no point saved): {row}")
    return row


def native_check(w, data, low, high):
    """Independent intensity NLL, including a continuous fixed-NP domain check.

    At fixed NPs each intensity is quadratic in sqrt(mu). Bounds alone do
    not guarantee positivity. Individual signed process terms are retained.
    """
    processes = ("S", "SBI", "B", "EWB", "EWSBI", "EWSBI10",
                 "qqZZ_0", "qqZZ_1", "qqZZ_2", "ttV")
    shapes = [w.function("prod_process_" + p) for p in processes]
    factors = [w.function("nfFactor_" + p) for p in processes]
    cr = [w.pdf(f"pdf_CR{c}") for c in range(3)]
    mu, original = w.var("mu"), w.var("mu").getVal()
    def rates():
        return np.array([f.getVal() for f in factors] +
                        [sum(f.getVal() for f in p.coefList()) for p in cr])
    current = rates()
    basis = []
    for value in (0., 1., 4.):
        mu.setVal(value)
        basis.append(rates())
    mu.setVal(original)
    basis = np.array(basis)
    variables, intensities, triples, weights = w.allVars(), [], [], []
    for i in range(data.numEntries()):
        row = data.get(i)
        channel = str(row.getCatLabel("cat"))
        if channel == "SR":
            variables.assignValueOnly(row)
            shape = np.array([f.getVal() for f in shapes])
            intensities.append(current[:10] @ shape)
            triples.append(basis[:, :10] @ shape)
        else:
            c = int(channel[-1])
            intensities.append(current[10+c])
            triples.append(basis[:, 10+c])
        weights.append(data.weight())
    intensities = np.array(intensities)
    if not np.all(np.isfinite(intensities)) or intensities.min() <= 0:
        raise RuntimeError("Nonpositive fitted intensity; no clipping is allowed")
    v0, v1, v2 = np.array(triples).T
    a, b = (v2 - 2*v1 + v0)/2, (4*v1 - 3*v0 - v2)/2
    left, right = math.sqrt(low), math.sqrt(high)
    minimum = np.minimum((a*left+b)*left+v0, (a*right+b)*right+v0)
    vertex = np.divide(-b, 2*a, out=np.full_like(a, left), where=a > 0)
    inside = (a > 0) & (vertex > left) & (vertex < right)
    minimum[inside] = np.minimum(minimum[inside], ((a*vertex+b)*vertex+v0)[inside])
    explicit = current.sum() - np.dot(weights, np.log(intensities))
    for alpha in w.allVars():
        name = str(alpha.GetName())
        if name.startswith("alpha_"):
            explicit += 0.5*(alpha.getVal() - w.var("nom_" + name[6:]).getVal())**2
    return dict(native_nll=float(explicit), minimum_intensity=float(intensities.min()),
                minimum_over_fixed_np_interval=float(minimum.min()))


def scan(args):
    ROOT.gROOT.SetBatch(True)
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
    if ROOT.gROOT.GetVersionInt() < 63404 or ROOT.gSystem.Load(str(args.library.resolve())) < 0:
        raise RuntimeError("Load the normalizer library built with ROOT >= 6.34.04 first")
    grid = sorted(set(args.grid))
    if not grid or any(not math.isfinite(x) or x < 0 or x > 2 for x in grid):
        raise ValueError("This demonstration supports finite mu values in [0, 2]")
    args.output.mkdir(parents=True, exist_ok=True)
    identity = dict(input_sha256=sha256(args.input), root_version=ROOT.gROOT.GetVersion(),
                    driver_sha256=sha256(__file__),
                    class_sha256=sha256(Path(__file__).parent / "normalizer/RooSameSampleNormalizer.cxx"),
                    class_header_sha256=sha256(Path(__file__).parent / "normalizer/RooSameSampleNormalizer.h"))
    manifest = json.loads(Path(__file__).with_name("workspaces.json").read_text())
    spec = next((s for s in manifest if s["sha256"] == identity["input_sha256"]), None)
    if spec is None:
        raise ValueError("Input does not match one of the three preserved workspaces")
    source = ROOT.TFile.Open(str(args.input), "READ")
    if not source or source.IsZombie():
        raise RuntimeError(f"Cannot open {args.input}")
    w = source.Get("combined")
    mc, data = w.obj("ModelConfig_4l_offshell"), w.data("obsData_4l_offshell")
    if not w.getSnapshot("asimov_generation"):
        raise RuntimeError("Missing Asimov generation snapshot")
    w.loadSnapshot("asimov_generation")
    pdf, mu = mc.GetPdf(), w.var("mu")
    parameters = pdf.getParameters(data)
    floating = [v for v in parameters if v.InheritsFrom("RooRealVar") and not v.isConstant()]
    expected = {v.GetName() for v in mc.GetNuisanceParameters()} | {v.GetName() for v in mc.GetParametersOfInterest()}
    if ({v.GetName() for v in floating} != expected or
        mc.GetNuisanceParameters().getSize() != spec["constrained_nps"] or
        not {"mu", "mu_qqZZ", "mu_qqZZ_1", "mu_qqZZ_2"}.issubset(expected)):
        raise RuntimeError("Unexpected floating parameters (preserve the pruning constants)")
    # Weighted extended Asimov likelihood, with the workspace's fixed globals.
    # Only mu is fixed during profiles; all three qqZZ factors stay free.
    nll = pdf.createNLL(data, ROOT.RooFit.Extended(True),
        ROOT.RooFit.Constrain(mc.GetNuisanceParameters()),
        ROOT.RooFit.GlobalObservables(mc.GetGlobalObservables()),
        ROOT.RooFit.GlobalObservablesSource("model"),
        ROOT.RooFit.EvalBackend("cpu"), ROOT.RooFit.Offset("initial"))
    nll.getVal()  # Initialize the same offset before restoring any checkpoint.
    generating = state(floating)
    domain = native_check(w, data, min(grid), max(grid))
    if domain["minimum_over_fixed_np_interval"] <= 0:
        raise RuntimeError("Generating nuisance point is not physical throughout this interval")
    checkpoint = args.output / "scan.json"
    if checkpoint.exists():
        saved = json.loads(checkpoint.read_text())
        if saved["identity"] != identity:
            raise RuntimeError("Checkpoint uses another input/ROOT/source; choose a new output folder")
    else:
        print(f"{spec['label']}: fitting {len(floating)} parameters from mu=1.4", flush=True)
        mu.setVal(1.4)  # Avoid the exactly stationary Asimov starting point.
        best = minimize(nll, 0.0001)
        best.update(mu_hat=mu.getVal(), parameters=state(floating))
        best.update(native_check(w, data, best["mu_hat"], best["mu_hat"]))
        saved = dict(identity=identity, label=spec["label"], profiled_nps=len(floating)-1,
                     generation_domain=domain, bestfit=best, points=[])
        save(checkpoint, saved)
    best = saved["bestfit"]
    restore(floating, best["parameters"])
    if not checked(best) or abs(nll.getVal() - best["nll"]) > 1e-6:
        raise RuntimeError("Checkpoint does not restore the global likelihood")
    print(f"{spec['label']}: mu_hat={best['mu_hat']:.9f}; profiling {len(floating)-1} NPs", flush=True)
    done = {round(r["mu"], 10): r for r in saved["points"]}
    branches = [sorted((x for x in grid if x < best["mu_hat"]), reverse=True),
                sorted(x for x in grid if x >= best["mu_hat"])]
    mu.setConstant(True)
    for branch in branches:
        restore(floating, best["parameters"])
        minim = make_minimizer(nll, 0.01)
        for value in branch:
            if round(value, 10) in done:
                restore(floating, done[round(value, 10)]["parameters"])
                continue
            mu.setVal(value)
            warm = state(floating)
            for fraction in (1., .5, .25, .125, 0.):
                for v in floating:
                    if v.GetName() != "mu":
                        name = v.GetName()
                        v.setVal(generating[name][0] + fraction*(warm[name][0]-generating[name][0]))
                if math.isfinite(nll.getVal()):
                    break
            else:
                raise RuntimeError(f"No physical start at mu={value}")
            if fraction != 1.:
                for v in floating:
                    if v.GetName() != "mu":
                        v.setError(generating[v.GetName()][1])
                minim = make_minimizer(nll, 0.01)
            print(f"Profiling mu={value:.3f}", flush=True)
            row = minimize(nll, 0.01, minim)
            row.update(mu=value, q=2*(row["nll"]-best["nll"]), parameters=state(floating))
            row.update(native_check(w, data, value, value))
            row["delta_nll_check"] = (row["native_nll"]-best["native_nll"]) - (row["nll"]-best["nll"])
            if row["q"] < -2e-4 or abs(row["delta_nll_check"]) > 1e-6:
                raise RuntimeError(f"Global minimum or independent NLL check failed at mu={value}")
            saved["points"].append(row)
            saved["points"].sort(key=lambda r: r["mu"])
            save(checkpoint, saved)
            print(f"mu={value:.3f}, q={row['q']:.8f}, EDM={row['edm']:.2g}, {row['seconds']:.1f} s", flush=True)
        del minim
    source.Close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--library", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--grid", type=float, nargs="+", required=True)
    scan(parser.parse_args())
