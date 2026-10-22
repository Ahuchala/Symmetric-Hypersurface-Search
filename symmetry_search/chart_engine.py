"""Reusable exact chart computation; no process management or checkpoints."""
import itertools
import time
from sage.all import QQ, GF, is_prime
from .models import coefficients, chart_polynomial


def coefficient_reduction(model,prime):
    K,terms=coefficients(model)
    rational=all(c in QQ for _,c in terms)
    if prime==0:
        return (QQ,QQ,{'characteristic':0}) if rational else (K,K,{'characteristic':0})
    if not is_prime(prime):raise ValueError('--prime must be 0 or prime')
    F=GF(prime)
    if rational:
        for _,c in terms:
            if QQ(c).denominator()%prime==0:raise ValueError('Prime divides a coefficient denominator')
        return F,lambda c:F(QQ(c)),{'characteristic':prime,'embedding':'rational'}
    if model['field_conductor']%prime==0:
        raise ValueError('Prime divides the cyclotomic conductor; choose another prime')
    roots=K.defining_polynomial().change_ring(F).roots(multiplicities=False)
    if not roots:raise ValueError('Coefficient field does not embed in F_p; choose a splitting prime')
    root=min(roots,key=int)
    for _,c in terms:
        if any(a.denominator()%prime==0 for a in c.list()):raise ValueError('Prime divides a coefficient denominator')
    return F,lambda c:c.polynomial().change_ring(F)(root),{'characteristic':prime,'cyclotomic_root':int(root)}


class ChartChecker:
    def __init__(self,model,prime,*,algorithm='auto',method='jacobian',sparse=True):
        if algorithm not in ('auto','std','slimgb'):raise ValueError('algorithm must be auto, std, or slimgb')
        if method not in ('jacobian','contraction'):raise ValueError('method must be jacobian or contraction')
        if method=='contraction' and (model['k'],model['d'])!=(3,1):
            raise ValueError('Contraction requires k=3 and d=1')
        self.model=model;self.algorithm=algorithm;self.method=method;self.sparse=sparse
        self.K,convert,self.embedding=coefficient_reduction(model,prime)
        self.terms=[(mon,convert(c)) for mon,c in coefficients(model)[1]]

    def equations(self,pivot):
        R,X,p,f=chart_polynomial(self.model,pivot,self.K,self.K,
                                sparse=self.sparse,terms=self.terms)
        if self.method=='jacobian':equations=[f]+[f.derivative(x) for x in R.gens()]
        else:
            # sigma(u_i,u_j,e_l)=0 for 3 pairs and all n coordinate vectors.
            # Equivalent to f=df=0 on the Grassmannian, in every characteristic.
            equations=[]
            for i,j in itertools.combinations(range(3),2):
                q=[R.zero() for _ in range(self.model['n'])]
                for mon,c in self.terms:
                    a,b,d=self.model['plucker_indices'][mon[0]]
                    q[a]+=c*(X[i,b]*X[j,d]-X[i,d]*X[j,b])
                    q[b]-=c*(X[i,a]*X[j,d]-X[i,d]*X[j,a])
                    q[d]+=c*(X[i,a]*X[j,b]-X[i,b]*X[j,a])
                equations.extend(q)
        return R,f,equations

    def check(self,pivot):
        start=time.perf_counter();R,f,equations=self.equations(pivot)
        built=time.perf_counter()
        result={'pivot':list(pivot),'embedding':self.embedding,'method':self.method,
                'algorithm':self.algorithm,'sparse_minors':self.sparse}
        if not f:
            result.update(status='invalid_reduction',reason='The whole section reduces to zero')
        else:
            I=R.ideal(equations)
            if self.algorithm=='auto':smooth=R.one() in I
            else:smooth=R.one() in I.groebner_basis(algorithm='libsingular:'+self.algorithm)
            result.update(status='smooth_chart' if smooth else 'singular_chart',
                          chart_polynomial_terms=len(f.dict()),singular_ideal_is_unit=bool(smooth))
        end=time.perf_counter()
        result['timings']={'equations_seconds':built-start,'ideal_seconds':end-built,'calculation_seconds':end-start}
        return result


def check_chart(model,pivot,prime,*,algorithm='auto',method='jacobian',sparse=True):
    return ChartChecker(model,prime,algorithm=algorithm,method=method,sparse=sparse).check(pivot)

def coordinate_singularity(checker):
    """Exact singular coordinate-plane witness, without a Groebner basis.

    On the p_I=1 chart, tangent Plucker coordinates are precisely p_J with
    |I intersection J|=k-1. Their signs do not affect vanishing. We require a
    nonzero standard-monomial section so zero reductions remain invalid,
    not misclassified as singular hypersurfaces.
    """
    model=checker.model;columns=[tuple(I) for I in model['plucker_indices']]
    lookup={I:j for j,I in enumerate(columns)};d=model['d'];k=model['k']
    polynomial={}
    for mon,c in checker.terms:
        mon=tuple(sorted(mon))
        if len(mon)!=d:return None
        if any(not all(a<=b for a,b in zip(columns[i],columns[j])) for i,j in zip(mon,mon[1:])):
            return None # conservative fallback for nonstandard representatives
        polynomial[mon]=polynomial.get(mon,checker.K.zero())+c
    polynomial={mon:c for mon,c in polynomial.items() if c}
    if not polynomial:return None
    for i,I in enumerate(columns):
        if polynomial.get((i,)*d):continue
        outside=set(range(model['n']))-set(I)
        neighbors=(lookup[tuple(sorted((set(I)-{a})|{b}))] for a in I for b in outside)
        if any(polynomial.get(tuple(sorted((i,)*(d-1)+(j,)))) for j in neighbors):continue
        return {'status':'singular_chart','pivot':list(I),'embedding':checker.embedding,
                'method':'coordinate_point','algorithm':'exact_evaluation',
                'singular_ideal_is_unit':False,
                'singular_point':{'nonzero_plucker_coordinate':list(I),'value':'1','chart_coordinates':'all zero'},
                'verification':'Nonzero standard-monomial section; coefficients of p_I^d and all tangent monomials p_I^(d-1)*p_J vanish exactly'}
    return None
