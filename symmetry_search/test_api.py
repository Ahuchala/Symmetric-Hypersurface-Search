"""Public interactive interface tests, run inside Sage."""
import contextlib
import io
import tempfile
import unittest
from pathlib import Path
from sage.all import QQ
from symmetry_search import search,load_search,load_hypersurface,CheckResult,Hypersurface
from .test_workflow import plain_model


class InteractiveAPI(unittest.TestCase):
    def test_in_memory_workflow_and_polynomials(self):
        output=io.StringIO()
        with contextlib.redirect_stdout(output):
            results=search(1,2,3,group='small:6:1',simple=False)
            X=results[1].construct(basis_tries=0)
            smooth=X.check_smoothness(seconds=20)
        self.assertEqual(output.getvalue(),'')
        self.assertTrue(results.complete)
        self.assertEqual(len(results),2)
        self.assertEqual(X.f.parent(),X.ring)
        self.assertEqual(X.f,X.p(0)**2*X.p(1)-X.p(0)*X.p(1)**2)
        self.assertEqual(len(X.invariant_basis),1)
        self.assertEqual(X.generators[0].nrows(),2)
        self.assertIs(smooth.is_smooth,True)
        self.assertIs(X.check_faithfulness().is_faithful,True)
        self.assertIs(results[results[1].id],results[1])
        with tempfile.TemporaryDirectory() as tmp:
            saved=Path(tmp)/'equation.json'
            X.save(saved)
            Y=load_hypersurface(saved)
            self.assertEqual(X.f,Y.f)
            # All charts were copied with the model; no new work is allowed here.
            self.assertIs(Y.check_smoothness(max_charts=0).is_smooth,True)
            results.save(Path(tmp)/'search.json')
            self.assertEqual(load_search(Path(tmp)/'search.json')[1].id,results[1].id)
        from symmetry_search import search as still_search
        self.assertTrue(callable(still_search))

    def test_partial_resume_and_no_ambiguous_truthiness(self):
        X=Hypersurface(plain_model(1,3,2,[((0,0),1),((1,1),1),((2,2),1)]))
        result=X.check_smoothness(prime=3,max_charts=1,seconds=20)
        self.assertIsNone(result.is_smooth)
        with self.assertRaises(TypeError):bool(result)
        self.assertIs(X.check_smoothness(prime=3,seconds=20).is_smooth,True)
        bad=X.check_smoothness(prime=2,seconds=20)
        self.assertEqual(bad.status,'singular_special_fibre')
        self.assertIsNone(bad.is_smooth)

    def test_coordinates_and_bounds(self):
        X=Hypersurface(plain_model(2,4,1,[((0,),1),((5,),1)]))
        self.assertEqual(X.p(1,0),-X.p(0,1))
        self.assertEqual(X.p(1,1),0)
        self.assertEqual(X.p((0,1)),X.p(0,1))
        with self.assertRaises(ValueError):X.p(4,0)
        with self.assertRaises(ValueError):search(1,2,3,simple='all')
        with self.assertRaises(ValueError):search(1,2,3,group=[])
        with self.assertRaises(ValueError):search(1,2,3,max_groups=0)
        with self.assertRaises(ValueError):X.check_smoothness(seconds=0)

    def test_incomplete_search_is_visible_and_data_is_copied(self):
        r=search(1,2,3,group=['small:6:1','small:2:1'],simple=None,max_groups=1)
        self.assertFalse(r.complete)
        self.assertIn('INCOMPLETE',repr(r))
        copy=r[0].data;copy['id']='modified'
        self.assertNotEqual(r[0].id,'modified')

    def test_stale_saved_equation_refused(self):
        X=Hypersurface(plain_model(1,2,2,[((0,0),1),((1,1),1)]))
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'equation.json';X.save(path)
            import json
            data=json.loads(path.read_text());data['terms'][0]['coefficient']=['2']
            path.write_text(json.dumps(data))
            with self.assertRaises(ValueError):X.check_smoothness(max_charts=0)

if __name__=='__main__':unittest.main()
