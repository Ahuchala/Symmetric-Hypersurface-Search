"""Independent small geometry oracles and orchestration failure regressions."""
import contextlib
from copy import deepcopy
import io
import itertools
import json
from pathlib import Path
import subprocess
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from sage.all import GF, QQ
from .benchmark import model, cases, compare_reports, run_benchmarks
from .common import write_json
from .models import chart_polynomial
from .smoothness import check_chart, coefficient_reduction, smooth


class GeometryRegression(unittest.TestCase):
    def test_singularity_only_at_boundary(self):
        equation=model(2,4,1,[((0,),1)])
        statuses={I:check_chart(equation,I,0)['status'] for I in itertools.combinations(range(4),2)}
        self.assertEqual(statuses[(0,1)],'smooth_chart') # empty hypersurface chart
        self.assertEqual([I for I,s in statuses.items() if s=='singular_chart'],[(2,3)])
        R,X,p,f=chart_polynomial(equation,(2,3),QQ,QQ)
        origin=[0]*R.ngens()
        self.assertEqual(f(*origin),0)
        self.assertTrue(all(f.derivative(x)(*origin)==0 for x in R.gens()))

    def test_finite_field_point_oracle(self):
        # Enumerate actual singular points independently of ideal membership.
        # For these rank-4/rank-2 alternating forms, the singular locus is
        # empty/a rational point, so rational enumeration detects it exactly.
        for terms in [[((0,),1),((5,),1)],[((0,),1)]]:
            equation=model(2,4,1,terms)
            for pivot in itertools.combinations(range(4),2):
                F=GF(3);R,X,p,f=chart_polynomial(equation,pivot,F,F)
                equations=[f]+[f.derivative(x) for x in R.gens()]
                witness=any(all(g(*point)==0 for g in equations)
                            for point in itertools.product(F,repeat=R.ngens()))
                self.assertEqual(check_chart(equation,pivot,3)['status'],
                                 'singular_chart' if witness else 'smooth_chart')

    def test_nonreduced_is_singular_not_reduced_support(self):
        reduced=model(1,2,1,[((0,),1)])
        doubled=model(1,2,2,[((0,0),1)])
        self.assertEqual(check_chart(reduced,(1,),0)['status'],'smooth_chart')
        self.assertEqual(check_chart(doubled,(1,),0)['status'],'singular_chart')

    def test_nonzero_scaling_preserves_characteristic_zero(self):
        original=model(2,4,1,[((0,),1),((5,),1)])
        scaled=model(2,4,1,[((0,),'-7/3'),((5,),'-7/3')])
        for pivot in itertools.combinations(range(4),2):
            self.assertEqual(check_chart(original,pivot,0)['status'],check_chart(scaled,pivot,0)['status'])
        with self.assertRaises(ValueError):coefficient_reduction(scaled,3)
        self.assertEqual(check_chart(scaled,(0,1),7)['status'],'invalid_reduction')

    def test_cyclotomic_embedding_and_invalid_primes(self):
        equation=model(1,2,1,[((0,),1),((1,),1)])
        equation['field_conductor']=3
        equation['terms'][1]['coefficient']=['0','1']
        self.assertEqual(check_chart(equation,(0,),0)['status'],'smooth_chart')
        self.assertEqual(check_chart(equation,(0,),7)['status'],'smooth_chart')
        for prime in [2,3,4]:
            with self.assertRaises(ValueError):coefficient_reduction(equation,prime)


class WorkerRegression(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.path=Path(self.tmp.name)/'model.json'
        write_json(self.path,model(1,2,1,[((0,),1)]))
        self.args=SimpleNamespace(model=str(self.path),prime=0,seconds=10,max_charts=2,
             output=str(Path(self.tmp.name)/'report.json'),restart=False,verbose=False)

    def run_check(self):
        with contextlib.redirect_stdout(io.StringIO()):return smooth(self.args)

    def test_timeout_is_inconclusive_and_retried(self):
        with patch('symmetry_search.smoothness.subprocess.run',side_effect=subprocess.TimeoutExpired('worker',10)):
            result=self.run_check()
        self.assertEqual(result['status'],'incomplete')
        self.assertTrue(all(c['status']=='timeout' for c in result['charts'].values()))
        self.assertEqual(self.run_check()['status'],'smooth_characteristic_zero')

    def test_worker_error_is_not_smooth(self):
        with patch('symmetry_search.smoothness.subprocess.run',return_value=SimpleNamespace(returncode=1,stderr='test failure')):
            result=self.run_check()
        self.assertEqual(result['status'],'incomplete')
        self.assertTrue(all(c['status']=='error' for c in result['charts'].values()))

    def test_completed_checkpoint_does_no_work(self):
        self.assertEqual(self.run_check()['status'],'smooth_characteristic_zero')
        with patch('symmetry_search.smoothness.subprocess.run',side_effect=AssertionError('recomputed cached chart')):
            self.assertEqual(self.run_check()['status'],'smooth_characteristic_zero')

    def test_prime_mismatch_rejected(self):
        self.run_check();self.args.prime=3
        with self.assertRaises(ValueError):self.run_check()

    def test_zero_budget_is_incomplete(self):
        self.args.max_charts=0
        result=self.run_check()
        self.assertEqual(result['status'],'incomplete')
        self.assertEqual(result['charts'],{})


class BenchmarkRegression(unittest.TestCase):
    def test_benchmark_failures_cannot_count_as_speedups(self):
        fixture={'id':'tiny','model':model(1,2,1,[((0,),1)]),
                 'prime':0,'pivot':[0],'expected':'smooth_chart'}
        with patch('symmetry_search.benchmark.platform.platform',return_value='test'), \
             patch('symmetry_search.benchmark.platform.processor',return_value='test'), \
             patch('symmetry_search.benchmark.cases',return_value=[fixture]), \
             patch('symmetry_search.benchmark.subprocess.check_output',return_value='test-revision'), \
             patch('symmetry_search.benchmark.subprocess.run',side_effect=subprocess.TimeoutExpired('worker',1)):
            result=run_benchmarks(repeats=1,seconds=1)
        self.assertFalse(result['passed'])
        self.assertEqual(result['cases'][0]['samples'][0]['status'],'timeout')
        self.assertNotIn('median_calculation_seconds',result['cases'][0])
        with self.assertRaises(ValueError):compare_reports(result,result)

    def test_baseline_cannot_be_overwritten(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'baseline.json';path.write_text('preserve me')
            with self.assertRaises(ValueError):run_benchmarks(output=path,baseline=path)
            self.assertEqual(path.read_text(),'preserve me')

    def test_full_has_every_a7_chart(self):
        a7=[c for c in cases('full') if c['id'].startswith('a7-')]
        self.assertEqual({tuple(c['pivot']) for c in a7},set(itertools.combinations(range(10),3)))
        self.assertEqual(len(a7),120)

    def test_comparison_rejects_changed_contract_and_failures(self):
        report={'schema':1,'suite_fingerprint':'fixed','settings':{'repeats':3},'passed':True,
                'environment':dict(sage='10',python='3',system='test',machine='test'),
                'cases':[dict(id='one',median_calculation_seconds=2,median_wall_seconds=4)]}
        faster=deepcopy(report);faster['cases'][0]['median_calculation_seconds']=1
        self.assertEqual(compare_reports(faster,report)[0]['calculation_speedup'],2)
        for key,value in [('passed',False),('suite_fingerprint','changed'),('settings',{})]:
            changed=deepcopy(report);changed[key]=value
            with self.assertRaises(ValueError):compare_reports(changed,report)
        changed=deepcopy(report);changed['environment']['sage']='different'
        with self.assertRaises(ValueError):compare_reports(changed,report)

if __name__=='__main__':unittest.main()
