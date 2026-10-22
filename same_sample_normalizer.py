"""Load the persistent RooFit normalizer with a source-addressed ACLiC cache."""
import hashlib
import fcntl
import os
from pathlib import Path

import ROOT


def load_normalizer():
    klass = ROOT.TClass.GetClass("RooSameSampleNormalizer", False, True)
    if klass and klass.HasDictionary():
        return
    source = Path(__file__).resolve().parent
    digest = hashlib.sha256(b"".join((source / name).read_bytes() for name in
                            ("RooSameSampleNormalizer.h", "RooSameSampleNormalizer.cxx"))).hexdigest()[:16]
    # Set this to a persistent /work directory for batch jobs. Never compile in
    # an input workspace directory or modify the installed ROOT distribution.
    build = Path(os.environ.get("OFFSHELL_ROOT_BUILD", "/tmp/offshell_root_build")) / digest
    build.mkdir(parents=True, exist_ok=True)
    library = build / "RooSameSampleNormalizer_cxx.so"
    with (build / "compile.lock").open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if library.exists():
            if ROOT.gSystem.Load(str(library)) < 0:
                raise RuntimeError(f"Cannot load {library}")
        else:
            ROOT.gSystem.SetBuildDir(str(build), True)
            if not ROOT.gSystem.CompileMacro(str(source / "RooSameSampleNormalizer.cxx"), "kO"):
                raise RuntimeError("Cannot compile RooSameSampleNormalizer")
    if not hasattr(ROOT, "RooSameSampleNormalizer"):
        raise RuntimeError("Missing RooSameSampleNormalizer dictionary")
