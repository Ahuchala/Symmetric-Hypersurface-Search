"""Run: sage -python -m unittest symmetry_search.test_workflow -v."""
import itertools
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from sage.all import QQ, GF, libgap, vector, identity_matrix, matrix, prod
from .search import initialize,candidates_for_table,explicit_table,constituents
from .sections import standard_monomials,section_dimension,invariant_space,verify_section,diagonal_simplification
from .models import chart_polynomial,make_model
from .representations import realize
from .faithfulness import check_faithfulness
from .smoothness import check_chart,smooth
from .common import encode,write_json,LimitExceeded


def plain_model(k,n,d,terms):
    cols=list(itertools.combinations(range(n),k))
    return {'schema':1,'k':k,'n':n,'d':d,'field_conductor':1,'plucker_indices':[list(c) for c in cols],
            'terms':[{'monomial':list(mon),'coefficient':encode(QQ(c))} for mon,c in terms]}

class Characters(unittest.TestCase):
    @classmethod
    def setUpClass(cls):initialize()

    def test_rectangles_and_trivial_rep_dimensions(self):
        tbl=libgap.TrivialGroup().CharacterTable();one=tbl.Irr()[0]
        for n in range(2,6):
            for k in range(1,n):
                for d in range(1,4):
                    chi=libgap.GrassSectionCharacter(tbl,n*one,k,d)
                    self.assertEqual(int(chi[0]),section_dimension(k,n,d))
                    _,mons=standard_monomials(k,n,d,10000)
                    self.assertEqual(len(mons),int(chi[0]))

    def test_plucker_relation_removed(self):
        tbl=libgap.CharacterTable('A5');chi=next(c for c in tbl.Irr() if int(c[0])==4)
        ext=libgap.GrassSectionCharacter(tbl,chi,2,1)
        ambient=libgap.SymmetricParts(tbl,[ext],2)[0]
        sections=libgap.GrassSectionCharacter(tbl,chi,2,2)
        one=tbl.TrivialCharacter()
        self.assertEqual(int(ambient[0]),21)
        self.assertEqual(int(sections[0]),20)
        self.assertEqual(int(libgap.ScalarProduct(tbl,ambient,one)),2)
        self.assertEqual(int(libgap.ScalarProduct(tbl,sections,one)),1)

    def test_scalar_kernel_is_not_faithful_projectively(self):
        tbl=libgap.DihedralGroup(8).CharacterTable()
        chi=next(c for c in tbl.Irr() if int(c[0])==2)
        self.assertEqual(int(libgap.GrassProjectiveKernelSize(tbl,chi)),2)
        self.assertEqual(list(candidates_for_table(tbl,'test',1,2,4)),[])

    def test_reducible_optional_and_bound(self):
        self.assertEqual(list(constituents([1,1],2,False,100)),[])
        self.assertEqual(list(constituents([1,1],2,True,100)),[[1,1],[1,2],[2,2]])
        with self.assertRaises(LimitExceeded):list(constituents([1,1],2,True,1))

    def test_nontrivial_linear_character(self):
        tbl,spec=explicit_table('small:6:1')
        found=list(candidates_for_table(tbl,spec,1,2,3))
        self.assertEqual(len(found),2)
        self.assertTrue(any(len(set(c['linear_character_values']))>1 for c in found))

class Geometry(unittest.TestCase):
    def test_smooth_and_singular_hyperplane_sections(self):
        # p01+p23 is a nondegenerate alternating form; p01 is decomposable.
        good=plain_model(2,4,1,[((0,),1),((5,),1)])
        bad=plain_model(2,4,1,[((0,),1)])
        self.assertTrue(all(check_chart(good,I,0)['status']=='smooth_chart' for I in itertools.combinations(range(4),2)))
        self.assertTrue(any(check_chart(bad,I,0)['status']=='singular_chart' for I in itertools.combinations(range(4),2)))

    def test_bad_prime_is_not_charzero_singularity(self):
        model=plain_model(1,3,2,[((0,0),1),((1,1),1),((2,2),1)])
        self.assertEqual(check_chart(model,(0,),0)['status'],'smooth_chart')
        self.assertEqual(check_chart(model,(0,),2)['status'],'singular_chart')
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'model.json';write_json(path,model)
            args=SimpleNamespace(model=str(path),prime=2,output=str(Path(tmp)/'report.json'),restart=False,max_charts=10,seconds=20)
            result=smooth(args)
            self.assertEqual(result['status'],'singular_special_fibre')
            self.assertIn('inconclusive',result['characteristic_zero_conclusion'])

    def test_partial_coverage_and_resume(self):
        model=plain_model(1,3,2,[((0,0),1),((1,1),1),((2,2),1)])
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'model.json';write_json(path,model)
            args=SimpleNamespace(model=str(path),prime=3,output=str(Path(tmp)/'report.json'),restart=False,max_charts=1,seconds=20)
            self.assertEqual(smooth(args)['status'],'incomplete')
            args.max_charts=10
            result=smooth(args)
            self.assertEqual(result['status'],'smooth_characteristic_zero')
            model['terms'][0]['coefficient']=['2'];write_json(path,model)
            with self.assertRaises(ValueError):smooth(args)

    def test_vanishing_mod_prime_and_denominators(self):
        model=plain_model(1,2,1,[((0,),3)])
        self.assertEqual(check_chart(model,(0,),3)['status'],'invalid_reduction')
        model=plain_model(1,2,1,[((0,),QQ(1)/3)])
        with self.assertRaises(ValueError):check_chart(model,(0,),3)

    def test_ambient_faithful_but_hypersurface_not(self):
        # C2 acts faithfully on P1, but fixes this single point.
        model=plain_model(1,2,1,[((0,),1)])
        model['classes']=[{'matrix':[[['1'],['0']],[['0'],['-1']]],'size':1,'order':2}]
        self.assertEqual(check_faithfulness(model,2)['status'],'not_faithful')
        # Two exchanged points give a faithful C2 action.
        model['d']=2;model['terms']=[{'monomial':[0,0],'coefficient':['1']},{'monomial':[1,1],'coefficient':['-1']}]
        self.assertEqual(check_faithfulness(model,2)['status'],'faithful')

    def test_diagonal_coefficient_simplification(self):
        cols,mons=standard_monomials(3,10,1)
        signed={(0,1,7):1,(0,4,6):-1,(0,5,9):1,(1,2,8):1,(1,5,6):-1,
                (2,3,9):1,(2,6,7):-1,(3,4,5):1,(3,7,8):-1,(4,8,9):-1}
        v=vector(QQ,[signed.get(I,0)*2**sum(j>=5 for j in I) for I in cols])
        S,new,basis=diagonal_simplification(QQ,cols,mons,v,[v])
        self.assertEqual({abs(c) for c in new if c},{1})
        ratios={new[i]*prod(S[j,j] for j in I)/v[i] for i,I in enumerate(cols) if v[i]}
        self.assertEqual(len(ratios),1)

    def test_nontrivial_character_explicit_construction(self):
        initialize();tbl,spec=explicit_table('small:6:1')
        candidate=next(c for c in candidates_for_table(tbl,spec,1,2,3) if len(set(c['linear_character_values']))>1)
        rep=realize(candidate);K=rep['field'];cols,mons=standard_monomials(1,2,3)
        basis=invariant_space(K,rep['generators'],rep['generator_characters'],1,2,3,cols,mons)
        self.assertEqual(len(basis),candidate['multiplicity'])
        verify_section(K,rep['generators'],rep['generator_characters'],1,2,3,cols,mons,basis[0])
        model=make_model(candidate,rep,'test',K,identity_matrix(K,2),cols,mons,basis[0],basis,[])
        self.assertEqual(check_faithfulness(model,2)['status'],'faithful')
        self.assertTrue(all(check_chart(model,I,0)['status']=='smooth_chart' for I in [(0,),(1,)]))

if __name__=='__main__':unittest.main()
