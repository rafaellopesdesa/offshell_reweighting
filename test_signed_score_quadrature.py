"""Equivalent positive/signed quadratures must give the same source variance."""
import unittest
import numpy as np
from compute_observation_response import score_integral


class SignedScoreTest(unittest.TestCase):
    def test_exact_quadrature_for_linear_scores(self):
        # Both rules integrate 1, x, x^2 to 2, 0, 1. The signed rule is
        # an integration representation of the same positive measure.
        x=np.array([-1.,0.,1.]); w=np.array([.5,1.,.5])
        signed_x=np.array([-1.,0.,1.,2.]); signed_w=np.array([-.5,4.,-2.5,1.])
        for offset,slope in [(0.,1.),(2.,3.),(-1.,.5)]:
            expected=2*offset**2+slope**2
            positive,_=score_integral(w,offset+slope*x)
            signed,parts=score_integral(signed_w,offset+slope*signed_x)
            self.assertAlmostEqual(positive,expected)
            self.assertAlmostEqual(signed,expected)
            self.assertAlmostEqual(parts['positive']+parts['negative'],expected)
        self.assertNotEqual(score_integral(abs(signed_w),signed_x)[0],1.)
        self.assertNotEqual(score_integral(signed_w**2,signed_x)[0],1.)

    def test_negative_estimate_is_not_silently_clipped(self):
        value,parts=score_integral([1.,-.5],[0.,2.])
        self.assertEqual(value,-2.)
        self.assertEqual(parts['negative'],-2.)


if __name__=='__main__':unittest.main()
