#!/usr/bin/env python3
"""Executable two-event example for RooSameSampleNormalizer (ROOT >= 6.34.04)."""
import argparse
import json
import math
import ROOT


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--library',help='Optional standalone shared library; otherwise use the project ACLiC loader')
    args=p.parse_args()
    if args.library:
        assert ROOT.gSystem.Load(args.library)>=0
    else:
        from same_sample_normalizer import load_normalizer
        load_normalizer()
    alpha=ROOT.RooRealVar('example_alpha','shape nuisance',0.,-3.,3.)
    vector=ROOT.std.vector('double')
    # u=(.4,.6), nominal rhat=(1,1), hence measure=u*rhat.
    # Arrays are flattened nuisance-major: K*N entries for each anchor.
    z=ROOT.RooSameSampleNormalizer('example_Z','same-reference integral',
        ROOT.RooArgList(alpha),vector([.4,.6]),vector([.9,1.1]),vector([1.2,.8]))
    rows=[]
    for value,expected in [(0.,1.),(1.,.96),(-1.,1.02),(1.,.96)]:
        alpha.setVal(value)
        assert math.isclose(z.getVal(),expected,rel_tol=2e-13)
        assert math.isclose(z.getVal(),z.directValue(),rel_tol=2e-13)
        rows.append({'alpha':value,'Z':z.getVal()})
    normalized=[1.2/z.getVal(),.8/z.getVal()]
    assert math.isclose(sum(u*f for u,f in zip([.4,.6],normalized)),1.,rel_tol=2e-13)
    # A separate coefficient lambda=100 retains exactly 100 expected events.
    assert math.isclose(sum(u*100*f for u,f in zip([.4,.6],normalized)),100.,rel_tol=2e-13)
    alpha.setVal(.25)
    assert math.isclose(z.getVal(),z.directValue(),rel_tol=2e-13)
    print(json.dumps({'passed':True,'root_version':ROOT.gROOT.GetVersion(),
        'anchors':rows,'normalized_high_shapes':normalized,'interior_Z':z.getVal()},indent=2))


if __name__=='__main__':main()
