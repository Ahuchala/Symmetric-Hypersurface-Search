import unittest
from unittest.mock import patch
from .api import search, SearchResults, Candidate, Hypersurface, CheckResult
from .common import LimitExceeded
from .smooth_search import combination


def smooth(status):return CheckResult({'status':status},'smoothness')
def faithful(status='faithful'):return CheckResult({'status':status},'faithfulness')

class SmoothDiscovery(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        results=search(1,2,2,group='small:2:1',simple=None,reducible=True)
        cls.c=next(c for c in results if c.multiplicity==2)
        cls.base=cls.c.construct(faithfulness_charts=0)

    def results(self,count=1):
        return SearchResults({'complete':True,'candidates':[self.c.data for _ in range(count)]})

    def test_reuses_basis_and_skips_singular_section(self):
        r=self.results()
        with patch.object(Candidate,'construct',return_value=self.base) as build, \
             patch.object(Hypersurface,'check_faithfulness',return_value=faithful()), \
             patch.object(Hypersurface,'check_smoothness',side_effect=[smooth('singular_characteristic_zero'),smooth('smooth_characteristic_zero')]):
            X=r.find_smooth(seed=7)
        self.assertIsNotNone(X);self.assertEqual(build.call_count,1)
        self.assertNotEqual(X.f,self.base.f)
        self.assertEqual(r.last_smooth_search['status'],'found')
        self.assertEqual(len(r.last_smooth_search['attempts']),2)
        self.assertEqual(self.base.data['verification']['smoothness'],'unchecked')

    def test_skips_failed_candidate_and_inconclusive_check(self):
        r=self.results(3)
        with patch.object(Candidate,'construct',side_effect=[LimitExceeded('budget'),self.base,self.base]), \
             patch.object(Hypersurface,'check_faithfulness',return_value=faithful()), \
             patch.object(Hypersurface,'check_smoothness',side_effect=[smooth('incomplete'),smooth('smooth_characteristic_zero')]):
            self.assertIsNotNone(r.find_smooth(sections_per_candidate=1))
        self.assertEqual([a['candidate_index'] for a in r.last_smooth_search['attempts']],[0,1,2])

    def test_no_false_success_and_candidate_bound(self):
        for status in ['not_faithful','inconclusive']:
            r=self.results(2)
            with patch.object(Candidate,'construct',return_value=self.base) as build, \
                 patch.object(Hypersurface,'check_faithfulness',return_value=faithful(status)), \
                 patch.object(Hypersurface,'check_smoothness') as check:
                self.assertIsNone(r.find_smooth(max_candidates=1,sections_per_candidate=1))
                self.assertEqual(build.call_count,1);check.assert_not_called()
            self.assertEqual(r.last_smooth_search['status'],'no_certificate_within_bounds')
        r=self.results()
        with patch.object(Candidate,'construct',return_value=self.base), \
             patch.object(Hypersurface,'check_faithfulness',return_value=faithful()), \
             patch.object(Hypersurface,'check_smoothness',return_value=smooth('singular_special_fibre')):
            self.assertIsNone(r.find_smooth(sections_per_candidate=1))

    def test_one_dimensional_space_and_interrupt(self):
        r=search(1,2,3,group='small:6:1',simple=False)
        X=r[1].construct(basis_tries=0,faithfulness_charts=0)
        with patch.object(Candidate,'construct',return_value=X), \
             patch.object(Hypersurface,'check_faithfulness',return_value=faithful()), \
             patch.object(Hypersurface,'check_smoothness',return_value=smooth('incomplete')) as check:
            self.assertIsNone(r.find_smooth(max_candidates=1,sections_per_candidate=5))
            self.assertEqual(check.call_count,1)
        with patch.object(Candidate,'construct',side_effect=KeyboardInterrupt):
            with self.assertRaises(KeyboardInterrupt):r.find_smooth()
        self.assertEqual(r.last_smooth_search['status'],'interrupted')

    def test_real_smooth_faithful_result_and_checkpoint(self):
        r=search(1,2,3,group='small:6:1',simple=False)
        X=r.find_smooth(prime=113,seconds=20,construction_options={'basis_tries':0})
        self.assertIsNotNone(X)
        self.assertIs(X.last_smoothness.is_smooth,True)
        self.assertIs(X.last_faithfulness.is_faithful,True)
        self.assertIs(X.check_smoothness(prime=113,max_charts=0).is_smooth,True)
        # The sparse C2-invariant square is unsuitable; a combination works.
        r=self.results()
        X=r.find_smooth(prime=113,seconds=20,seed=0)
        self.assertIsNotNone(X)
        self.assertGreaterEqual(len(r.last_smooth_search['attempts']),2)
        self.assertIs(X.last_smoothness.is_smooth,True)
        self.assertIs(X.last_faithfulness.is_faithful,True)

    def test_combination_and_validation(self):
        X=combination(self.base,[2,3])
        self.assertEqual(X.f,self.base.invariant_basis[0]+3/2*self.base.invariant_basis[1])
        self.assertIsNone(X.last_smoothness)
        with self.assertRaises(ValueError):combination(self.base,[0,0])
        for kwargs in [{'sections_per_candidate':0},{'seconds':0},{'construction_options':{'coefficients':[1]}},{'algorithm':'bad'}]:
            with self.assertRaises(ValueError):self.results().find_smooth(**kwargs)

if __name__=='__main__':unittest.main()
