#!/usr/bin/env python3
"""Stage small, seeded reference subsets from read-only NumPy inputs on Unity."""

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import time

import numpy as np


PROCESSES = ("S", "SBI", "B", "EWB", "EWSBI", "EWSBI10", "qqZZ", "ttV")
SR_PROCESSES = tuple(p for p in PROCESSES if p != "qqZZ") + ("qqZZ_0", "qqZZ_1", "qqZZ_2")
YIELD_INDEX = {"S": 0, "SBI": 1, "B": 2, "qqZZ_0": 3, "qqZZ_1": 4,
               "qqZZ_2": 5, "qqZZ": 6, "ttV": 7, "EWB": 8, "EWSBI": 9, "EWSBI10": 10}
CR_ROWS = ("SBI", "S", "B", "qqZZ", "ttV", "EWB", "EWSBI", "EWSBI10")


def write_json(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def sha256(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def nuisance_name(name):
    aliases = {"ATLAS_LUMI": "alpha_LUMI", "gg_PDF_var": "gg_PDF",
               "EW_PDF_var": "VBF_PDF", "HOEW_syst": "HOEW",
               "H4l_Shower_UEPS_Sherpa_HM_QSF": "PS_qqZZ_QSF_Shape",
               "ggZZNLO_QCD_syst_shape": "QCD_ggZZk_Shape",
               "ggZZNLO_QCD_syst_norm": "QCD_ggZZk_Norm"}
    for jet in range(3):
        aliases[f"qq_PDF_var_{jet}jet"] = f"qq_PDF_{jet}Jet"
        aliases[f"HOQCD_{jet}jet"] = f"HOQCD_{jet}Jet"
        aliases[f"HOEW_QCD_syst_{jet}jet"] = f"HOEW_QCD_{jet}Jet"
        aliases[f"{jet}jet_H4l_Shower_UEPS_Sherpa_HM_CKKW"] = f"PS_qqZZ_CKKW_{jet}Jet_Shape"
    for variation in ("CKKW", "QSF", "CSSKIN"):
        for part in ("norm", "shape"):
            aliases[f"H4l_Shower_UEPS_Sherpa_ggHM_{variation}_{part}"] = f"PS_ggZZ_{variation}_{part.title()}"
    for kind in ("Composition", "Response"):
        for source, target in (("ggF", "gg"), ("EW", "VBF"), ("qqZZ", "qq")):
            aliases[f"JET_Flavor_{kind}_{source}"] = f"JET_Flavor_{kind}_{target}"
    return "alpha_ATLAS_" + aliases.get(name, name)


def stage(source, output, maximum, seed, workers):
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Stage the large systematic arrays inside a Slurm allocation")
    if not 1 <= workers <= int(os.environ.get("SLURM_CPUS_ON_NODE", "1")):
        raise ValueError("I/O workers must fit within the Slurm CPU allocation")
    output.mkdir(parents=True, exist_ok=True)
    if (output / "stage_manifest.json").exists():
        raise FileExistsError("Staging already completed; reuse it or select a new output directory")
    weight = np.load(source / "w_SM.npy", mmap_mode="r", allow_pickle=False)
    count = len(weight)
    if not 0 < maximum <= count:
        raise ValueError("Invalid subset size")
    indices = np.random.default_rng(seed).permutation(count)[:maximum]
    # Sort disk reads, then restore random-permutation order in the saved slices.
    order = np.argsort(indices)
    sorted_indices = indices[order]
    inverse = np.argsort(order)
    index_path = output / "event_rows.npy"
    if index_path.exists():
        np.testing.assert_array_equal(np.load(index_path, mmap_mode="r"), indices)
    else:
        np.save(index_path, indices)
    progress_path = output / "stage_progress.json"
    progress = json.loads(progress_path.read_text()) if progress_path.exists() else {}
    for key, value in (("seed", seed), ("maximum_subset", maximum), ("population", count)):
        if key in progress and progress[key] != value:
            raise ValueError(f"Staging checkpoint has a different {key}")
        progress[key] = value
    inputs, outputs = progress.setdefault("inputs", {}), progress.setdefault("outputs", {})
    def save(name, array):
        path = output / name
        if array.dtype.kind not in "US" and not np.isfinite(array).all():
            raise ValueError(f"Non-finite input in {name}")
        temporary = path.with_suffix(".partial.npy")
        np.save(temporary, array)
        temporary.replace(path)
        return {"sha256": sha256(path), "shape": list(array.shape), "dtype": str(array.dtype)}

    names = ["w_SM.npy", "njets.npy", "list_syst.npy", "yields_array_SR.npy", "yields_array_CR.npy", "hist_combined.npy"]
    names += [f"r_{p}vsRef.npy" for p in PROCESSES]
    names += [f"list_weight_{d}_{p}.npy" for p in SR_PROCESSES for d in ("down", "up")]
    names += [f"arr_ratio_{d}_{p}_CR.npy" for p in PROCESSES for d in ("down", "up")]
    names += [f"arr_ratio_{d}_{p}.npy" for p in PROCESSES for d in ("down", "up")]
    def extract(name):
        path = source / name
        stat = path.stat()
        started = time.monotonic()
        array = np.load(path, mmap_mode="r", allow_pickle=False)
        record = {"path": str(path.resolve()), "bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns,
                  "shape": list(array.shape), "dtype": str(array.dtype)}
        if name in inputs and name in outputs:
            for key in record:
                if record[key] != inputs[name][key]:
                    raise ValueError(f"Source changed since checkpoint: {name}")
            if sha256(output / name) != outputs[name]["sha256"]:
                raise ValueError(f"Staged slice changed since checkpoint: {name}")
            return name, inputs[name], outputs[name], 0.0, True
        if array.ndim == 2 and array.shape[1] == count:
            if array.shape[0] != 127:
                raise ValueError(f"Unexpected systematic dimensions: {name}")
            selected = np.empty((127, maximum), dtype=array.dtype)
            for row in range(127):
                selected[row] = array[row, sorted_indices][inverse]
        elif array.shape == (count,):
            selected = array[sorted_indices][inverse]
        else:
            selected = np.array(array)
            record["sha256"] = sha256(path)
        saved = save(name, selected)
        after = path.stat()
        if (after.st_size, after.st_mtime_ns) != (stat.st_size, stat.st_mtime_ns):
            raise RuntimeError(f"Input changed while staging: {path}")
        return name, record, saved, time.monotonic()-started, False
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = [executor.submit(extract, name) for name in names]
        for future in as_completed(futures):
            name, record, saved, seconds, reused = future.result()
            inputs[name], outputs[name] = record, saved
            write_json(progress_path, progress)
            action = "Reused" if reused else "Staged"
            print(f"{action} {name}: {saved['shape']}, {seconds:.1f} s", flush=True)
    outputs["event_rows.npy"] = {"sha256": sha256(output / "event_rows.npy"), "shape": [maximum], "dtype": str(indices.dtype)}
    write_json(output / "stage_manifest.json", {
        "created_utc": datetime.now(timezone.utc).isoformat(), "slurm_job_id": os.environ["SLURM_JOB_ID"],
        "seed": seed, "population": count, "maximum_subset": maximum,
        "selection": "numpy.default_rng(seed).permutation(population)[:maximum]; no replacement",
        "numpy_version": np.__version__, "source": str(source.resolve()), "inputs": inputs, "outputs": outputs,
        "io_workers": workers,
        "integrity": "All inputs opened with mmap_mode='r'; source sizes and mtimes unchanged during reads. Large arrays are not fully hashed; selected slices are hashed.",
    })


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--maximum", type=int, default=8192)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()
    stage(args.source, args.output, args.maximum, args.seed, args.workers)


if __name__ == "__main__":
    main()
