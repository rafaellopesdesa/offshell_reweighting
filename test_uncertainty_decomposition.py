"""A soluble Gaussian experiment checks source closure and impact semantics."""
import json
import hashlib
from pathlib import Path
import tempfile
import unittest

import numpy as np

from report_uncertainty_decomposition import analyze


class SourcePropagationTest(unittest.TestCase):
    def test_measured_nuisance_and_unconstrained_background(self):
        # x=mu+alpha+gamma (var4), b=gamma (var9), y=alpha (var4),
        # and an auxiliary measurement a=alpha (var1).
        # alpha_hat=.8*a+.2*y; mu_hat=x-b-alpha_hat.
        names=["alpha_ATLAS_alpha_LUMI","mu","mu_qqZZ"]
        designs=[np.array([1.,1.,1.]),np.array([0.,0.,1.]),np.array([1.,0.,0.])]
        matrices={k:np.outer(v,v)/variance for k,v,variance in zip(("SR","CR0","CR1"),designs,(4.,9.,4.))}
        G=sum(matrices.values())+np.diag([1.,0.,0.]);V=np.linalg.inv(G)
        h={"passed":True,"input_sha256":"toy","parameter_order":names,"covariance":V.tolist(),
           "parameters":[{"name":n,"value":0.} for n in names],"covariance_quality":3,"hesse_status":0,"seconds":0.}
        b={"passed":True,"input_sha256":"toy","parameter_order":names,"parameters":dict.fromkeys(names,0.),
           "auxiliaries":{names[0]:{"value":0.}},"primary_variance":{k:float(V[1]@v@V[:,1]) for k,v in matrices.items()},
           "step_relative_variance_differences":{"SR":0.}}
        with tempfile.TemporaryDirectory() as tmp:
            directory=Path(tmp)
            (directory/"covariance.json").write_text(json.dumps(h))
            covariance_hash=hashlib.sha256((directory/"covariance.json").read_bytes()).hexdigest()
            b["covariance_sha256"]=covariance_hash
            (directory/"observation_response.json").write_text(json.dumps(b))
            (directory/"curvature_validation.json").write_text(json.dumps({"passed":True,"covariance_sha256":covariance_hash,
                "checks":[{"direction":"profile_mu_tangent","richardson_curvature":1.}]}))
            r=analyze(directory)
        self.assertAlmostEqual(r["hessian_variance"],13.8)
        self.assertAlmostEqual(r["auxiliary_variance"],.64)
        self.assertAlmostEqual(r["primary_variance"],13.16)
        self.assertAlmostEqual(r["closure_relative"],0.)
        self.assertAlmostEqual(r["conditional_diagnostics"]["alpha_fixed_sigma"]**2,13.)
        self.assertAlmostEqual(r["conditional_diagnostics"]["qqzz_block_profiling_variance_penalty"],9.)
        alpha=next(x for x in r["impacts"] if x["constrained"])
        self.assertAlmostEqual(alpha["unit_auxiliary_response"],-.8)
        self.assertAlmostEqual(alpha["conditional_postfit_np_impact"],-np.sqrt(.8))
        qq=next(x for x in r["impacts"] if not x["constrained"])
        self.assertIsNone(qq["auxiliary_variance"])


if __name__=="__main__":unittest.main()
