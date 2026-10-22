"""Cross-check optimized algebra, symmetry certificates, and process lifecycle."""
import contextlib
import io
import itertools
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch
from sage.all import QQ
from .benchmark import model
from .chart_engine import ChartChecker, check_chart, coordinate_singularity
from .chart_symmetry import chart_orbits, verified_action
from .common import ROOT, write_json
from .models import chart_polynomial
from .smoothness import smooth
from .workers import PersistentWorker


class OptimizedAlgebra(unittest.TestCase):
    def test_coordinate_precheck_matches_actual_chart_derivatives(self):
        models=[model(2,4,1,[((0,),1)]),
                model(2,4,1,[((0,),1),((5,),1)]),
                json.loads((ROOT/'output/symmetry-examples/a5-gr25-cubic.json').read_text())]
        for m in models:
            checker=ChartChecker(m,0);witness=coordinate_singularity(checker)
            actual=[]
            for pivot in itertools.combinations(range(m['n']),m['k']):
                R,f,equations=checker.equations(pivot)
                if all(g(*([0]*R.ngens()))==0 for g in equations):actual.append(pivot)
            self.assertEqual(witness is not None,bool(actual))
            if witness:self.assertIn(tuple(witness['pivot']),actual)
        self.assertIn((2,3),actual)

    def test_precheck_skips_zero_and_nonstandard_sections(self):
        self.assertIsNone(coordinate_singularity(ChartChecker(model(1,2,1,[((0,),3)]),3)))
        # The Plucker relation is the zero section, despite nonzero coefficients.
        relation=model(2,4,2,[((0,5),1),((1,4),-1),((2,3),1)])
        self.assertIsNone(coordinate_singularity(ChartChecker(relation,0)))

    def test_a5_singularity_needs_no_worker(self):
        path=ROOT/'output/symmetry-examples/a5-gr25-cubic.json'
        with tempfile.TemporaryDirectory() as tmp:
            args=SimpleNamespace(model=str(path),prime=0,seconds=20,max_charts=10,
                output=str(Path(tmp)/'report.json'),restart=True,verbose=False,
                backend='persistent',workers=1,algorithm='slimgb',symmetry=True)
            with patch.object(PersistentWorker,'check',side_effect=AssertionError('unnecessary Groebner worker')),contextlib.redirect_stdout(io.StringIO()):
                result=smooth(args)
        self.assertEqual(result['status'],'singular_characteristic_zero')
        self.assertEqual(result['last_run']['worker_starts'],0)
        self.assertIn('singular_point',result['coordinate_singularity_precheck'])

    def test_sparse_minors_equal_dense_equation(self):
        m=json.loads((ROOT/'output/a7-equation.json').read_text())
        for pivot in [(0,1,2),(1,4,7),(7,8,9)]:
            dense=chart_polynomial(m,pivot,sparse=False)
            sparse=chart_polynomial(m,pivot,sparse=True)
            self.assertEqual(dense[3],sparse[3])
            self.assertEqual(sum(p is not None for p in sparse[2]),12)

    def test_contraction_ideal_equals_jacobian_ideal(self):
        # Includes proper (singular) ideals: equality cannot pass merely because
        # both sides happen to be unit ideals on a smooth example.
        for terms in [[((0,),1)],[((0,),1),((5,),1)]]:
            m=model(3,5,1,terms)
            for prime in [0,2,3,7]:
                J=ChartChecker(m,prime,method='jacobian')
                C=ChartChecker(m,prime,method='contraction')
                for pivot in itertools.combinations(range(5),3):
                    R,f,eqJ=J.equations(pivot);_,_,eqC=C.equations(pivot)
                    IJ=R.ideal(eqJ);IC=R.ideal(eqC)
                    self.assertTrue(all(g in IC for g in eqJ))
                    self.assertTrue(all(g in IJ for g in eqC))
                    self.assertTrue(all(g.total_degree()<=2 for g in eqC if g))

    def test_algorithm_agreement(self):
        fixtures=[model(2,4,1,[((0,),1),((5,),1)]),model(2,4,1,[((0,),1)])]
        for m in fixtures:
            for pivot in itertools.combinations(range(4),2):
                results=[check_chart(m,pivot,7,algorithm=a)['status'] for a in ['auto','std','slimgb']]
                self.assertEqual(len(set(results)),1)

    def test_reject_invalid_options(self):
        m=model(2,4,1,[((0,),1)])
        with self.assertRaises(ValueError):ChartChecker(m,0,method='contraction')
        with self.assertRaises(ValueError):ChartChecker(m,0,algorithm='unknown')

    def test_verified_symmetry_and_rejected_guess(self):
        m=model(2,4,1,[((0,),1),((5,),1)])
        orbits,proofs=chart_orbits(m,0,permutations=[[2,3,0,1]])
        self.assertLess(len(orbits),6)
        self.assertEqual(sum(map(len,orbits)),6)
        with self.assertRaises(ValueError):chart_orbits(m,0,permutations=[[1,0,2,3]])
        self.assertIsNone(verified_action(m,[0,1,2,3],[QQ(7)]*4,7))

    def test_a7_support_symmetry_covers_all_charts(self):
        m=json.loads((ROOT/'output/a7-equation.json').read_text())
        orbits,proofs=chart_orbits(m,113,True)
        self.assertEqual(len(orbits),42)
        self.assertEqual({I for orbit in orbits for I in orbit},set(itertools.combinations(range(10),3)))
        self.assertTrue(proofs)


class OptimizedWorkers(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.path=Path(self.tmp.name)/'equation.json'
        write_json(self.path,model(2,4,1,[((0,),1),((5,),1)]))

    def test_reuse_timeout_kill_and_replacement(self):
        with PersistentWorker(self.path,0) as worker:
            self.assertEqual(worker.check((0,1),20)['status'],'smooth_chart')
            process=worker.process
            self.assertEqual(worker.check((0,2),20)['status'],'smooth_chart')
            self.assertIs(worker.process,process);self.assertEqual(worker.starts,1)
            with patch.object(worker,'_line',side_effect=TimeoutError):
                self.assertEqual(worker.check((0,3),20)['status'],'timeout')
            self.assertIsNotNone(process.poll())
            self.assertEqual(worker.check((1,2),20)['status'],'smooth_chart')
            self.assertEqual(worker.starts,2)
        self.assertIsNone(worker.process)

    def test_real_timeout_leaves_no_live_process(self):
        with PersistentWorker(self.path,0) as worker:
            self.assertEqual(worker.check((0,1),0.000001)['status'],'timeout')
            self.assertIsNone(worker.process)

    def test_crash_recovery(self):
        with PersistentWorker(self.path,0) as worker:
            worker.check((0,1),20)
            worker.process.kill();worker.process.wait()
            self.assertEqual(worker.check((0,2),20)['status'],'error')
            self.assertEqual(worker.check((0,2),20)['status'],'smooth_chart')

    def test_parallel_symmetry_and_checkpoint(self):
        args=SimpleNamespace(model=str(self.path),prime=0,seconds=20,max_charts=100,
            output=str(Path(self.tmp.name)/'report.json'),restart=False,verbose=False,
            backend='persistent',workers=2,symmetry=True,algorithm='slimgb')
        with contextlib.redirect_stdout(io.StringIO()):result=smooth(args)
        self.assertEqual(result['status'],'smooth_characteristic_zero')
        self.assertEqual(len(result['charts']),6)
        self.assertTrue(any('inferred_from' in c for c in result['charts'].values()))
        self.assertLessEqual(result['last_run']['worker_starts'],2)
        with patch.object(PersistentWorker,'check',side_effect=AssertionError('recomputed')),contextlib.redirect_stdout(io.StringIO()):
            resumed=smooth(args)
        self.assertEqual(resumed['status'],'smooth_characteristic_zero')
        self.assertEqual(resumed['last_run']['charts_computed'],0)

if __name__=='__main__':unittest.main()
