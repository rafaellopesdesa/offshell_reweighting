"""Compile a shared RooFit library and its ROOT I/O dictionary."""
import argparse
import os
from pathlib import Path
import shlex
import subprocess

import ROOT


def build(directory):
    if ROOT.gROOT.GetVersionInt() < 63404:
        raise RuntimeError("ROOT >= 6.34.04 is required")
    directory = Path(directory).resolve()
    directory.mkdir(parents=True, exist_ok=True)
    source = Path(__file__).resolve().parent
    root = Path(str(ROOT.gROOT.GetRootSys()))
    include, lib = root / "include", root / "lib"
    config = subprocess.check_output([str(root / "bin/root-config"), "--cflags"], text=True)
    # Wheel root-config can contain build-time include paths. Use its active
    # installation's headers, retaining the matching C++ standard and flags.
    flags = [flag for flag in shlex.split(config) if not flag.startswith("-I")]
    library = directory / "libRooSameSampleNormalizer.so"
    dictionary = directory / "G__RooSameSampleNormalizer.cxx"
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = str(lib) + os.pathsep + env.get("LD_LIBRARY_PATH", "")
    subprocess.run([str(root / "bin/rootcling"), "-f", str(dictionary),
                    "-s", str(library), f"-I{include}", f"-I{source}",
                    str(source / "RooSameSampleNormalizer.h"),
                    str(source / "LinkDef.h")], check=True, env=env)
    subprocess.run([os.environ.get("CXX", "c++"), *flags, "-O2", "-fPIC", "-shared",
                    f"-I{include}", f"-I{source}",
                    str(source / "RooSameSampleNormalizer.cxx"), str(dictionary),
                    f"-L{lib}", f"-Wl,-rpath,{lib}",
                    "-lRooFit", "-lRooFitCore", "-lCore", "-lRIO", "-lHist",
                    "-lTree", "-lMathCore", "-o", str(library)], check=True, env=env)
    if ROOT.gSystem.Load(str(library)) < 0 or not ROOT.TClass.GetClass("RooSameSampleNormalizer").HasDictionary():
        raise RuntimeError("Cannot load the shared library and its I/O dictionary")
    print(f"ROOT {ROOT.gROOT.GetVersion()}: {library}", flush=True)
    return library


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("build/normalizer"))
    build(parser.parse_args().output)
