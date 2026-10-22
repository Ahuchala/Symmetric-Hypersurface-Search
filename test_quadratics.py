"""Exact cross-checks of Reynolds quadrics and generic Plucker straightening."""
import itertools
import unittest
from sage.all import QQ, ZZ, CyclotomicField, matrix, diagonal_matrix, vector, identity_matrix
from .common import LimitExceeded
from .sections import standard_monomials,invariant_space,verify_section
from .quadratics import QuadraticStraightener,invariant_quadrics,reynolds_gram
from .representations import realize
from . import search


class QuadraticConstruction(unittest.TestCase):
    def test_all_straightening_identities_in_several_grassmannians(self):
        for k,n in [(1,4),(2,4),(2,5),(3,5),(4,5)]:
            columns,monomials=standard_monomials(k,n,2,10000)
            S=QuadraticStraightener(k,n,columns,monomials,2000000)
            for pair in itertools.combinations_with_replacement(range(len(columns)),2):
                expression=S.reduce(pair);result={}
                for index,c in expression.items():
                    for e,a in S.polynomial(monomials[index]).items():result[e]=result.get(e,0)+c*a
                self.assertEqual({e:a for e,a in result.items() if a},S.polynomial(pair))

    def test_plucker_relation_reduces_to_zero(self):
        columns,monomials=standard_monomials(2,4,2)
        S=QuadraticStraightener(2,4,columns,monomials,10000)
        v=vector(QQ,len(monomials))
        for pair,c in [((0,5),1),((1,4),-1),((2,3),1)]:
            for j,a in S.reduce(pair).items():v[j]+=c*a
        self.assertFalse(v)

    def test_cyclotomic_gram_and_scalar_weights(self):
        K=CyclotomicField(3);z=K.gen()
        V=matrix(K,[[1,z,z+1],[z*z,2,z],[1-z,1,z*z]])
        for weights in [[K.one()]*3,[1,z,z*z]]:
            self.assertEqual(reynolds_gram(V,weights),V.transpose()*diagonal_matrix(K,weights)*V)

    def test_reynolds_quotient_agrees_with_expanded_a5(self):
        c=search(2,4,2,group='A5')[0].data;rep=realize(c)
        cols,mons,basis,evidence=invariant_quadrics(rep,2,4,c['multiplicity'])
        old=invariant_space(rep['field'],rep['generators'],rep['generator_characters'],2,4,2,cols,mons)
        self.assertEqual(matrix(rep['field'],basis).row_space(),matrix(rep['field'],old).row_space())
        self.assertEqual(len(basis),1) # ambient invariant quadrics have dimension 2
        self.assertEqual(evidence['evaluation_rank'],1)
        verify_section(rep['field'],rep['generators'],rep['generator_characters'],2,4,2,cols,mons,basis[0])

    def test_nontrivial_character_and_api(self):
        results=search(1,2,2,group='small:2:1',simple=None,reducible=True)
        c=next(c for c in results if len(set(c.data['linear_character_values']))>1)
        X=c.construct(faithfulness_charts=0)
        self.assertEqual(X.f,X.p(0)*X.p(1))
        self.assertEqual(X.data['quadratic_construction']['invariant_dimension'],1)
        self.assertEqual(c.construct(quadratic_method='expanded',faithfulness_charts=0).f,X.f)
        self.assertIs(X.check_faithfulness().is_faithful,False)
        with self.assertRaises(LimitExceeded):c.construct(max_sections=1)
        trivial=next(t for t in results if t.multiplicity==2)
        combination=trivial.construct(coefficients=[2,3],faithfulness_charts=0)
        self.assertEqual(combination.f,2*combination.p(0)**2+3*combination.p(1)**2)

    def test_number_field_representation(self):
        K=CyclotomicField(3);z=K.gen()
        rep={'field':K,'generators':[diagonal_matrix(K,[1,z,z*z])],
             'generator_characters':[K.one()],'group_order':3}
        cols,mons,basis,evidence=invariant_quadrics(rep,2,3,2)
        old=invariant_space(K,rep['generators'],rep['generator_characters'],2,3,2,cols,mons)
        self.assertEqual(matrix(K,basis).row_space(),matrix(K,old).row_space())
        self.assertEqual(evidence['evaluation_rank'],2)

    def test_nonrational_character_with_rational_matrices(self):
        K=CyclotomicField(3);z=K.gen()
        rep={'field':K,'generators':[matrix(K,[[0,-1],[1,-1]])],
             'generator_characters':[z],'group_order':3}
        cols,mons,basis,evidence=invariant_quadrics(rep,1,2,1)
        old=invariant_space(K,rep['generators'],rep['generator_characters'],1,2,2,cols,mons)
        self.assertEqual(matrix(K,basis).row_space(),matrix(K,old).row_space())
        self.assertEqual(evidence['coefficient_subfield_degree'],1)

    def test_denominators_in_group_orbits(self):
        C=matrix(QQ,[[1,2],[1,5]])
        M=C.inverse()*matrix(QQ,[[0,-1],[1,-1]])*C
        self.assertNotEqual(M.denominator(),1)
        rep={'field':QQ,'generators':[M],'generator_characters':[QQ.one()],'group_order':3}
        cols,mons,basis,_=invariant_quadrics(rep,1,2,1)
        old=invariant_space(QQ,rep['generators'],rep['generator_characters'],1,2,2,cols,mons)
        self.assertEqual(matrix(QQ,basis).row_space(),matrix(QQ,old).row_space())

    def test_resource_and_trial_bounds(self):
        rep={'field':QQ,'generators':[diagonal_matrix(QQ,[1,-1])],
             'generator_characters':[QQ.one()],'group_order':2}
        with self.assertRaises(LimitExceeded):invariant_quadrics(rep,1,2,2,max_entries=1)
        with self.assertRaises(LimitExceeded):invariant_quadrics(rep,1,2,2,max_trials=1)
        _,_,a,_=invariant_quadrics(rep,1,2,2,seed=0)
        _,_,b,_=invariant_quadrics(rep,1,2,2,seed=7)
        self.assertEqual(a,b)

if __name__=='__main__':unittest.main()
