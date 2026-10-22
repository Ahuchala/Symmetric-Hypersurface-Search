"""Exact invariant quadrics via Reynolds forms and Plucker straightening.

The computation is uniform in (k,n). It constructs symmetric forms on
Lambda^k(W), then certifies independence in the Grassmannian coordinate ring
by exact evaluations on decomposable vectors. No expanded invariance system
on S_(2^k)(W) is built. Character multiplicity certifies completeness.
"""
import itertools
import random
import time
from sage.all import QQ, ZZ, PolynomialRing, matrix, vector, identity_matrix, factorial, VectorSpace, block_matrix, lcm, prod
from .common import LimitExceeded, encode_matrix
from .sections import (section_dimension, standard_monomials, compound, minors,
                       normalize, choose_vector, diagonal_simplification, verify_section)


class QuadraticStraightener:
    """Subduction in a generic k by n matrix using diagonal leading terms.

    Each subtraction is an exact integer polynomial identity. Standard
    two-column tableaux have distinct diagonal leading monomials, giving
    a unique normal form without a Groebner computation in Plucker variables.
    """
    def __init__(self,k,n,columns,monomials,max_entries):
        self.k=k;self.n=n;self.columns=columns;self.limit=max_entries;self.entries=0
        self.column_index={I:j for j,I in enumerate(columns)}
        self.index={m:j for j,m in enumerate(monomials)}
        if int(factorial(k))**2>max_entries:
            raise LimitExceeded('Generic-minor straightening exceeds max_entries; increase that bound for this k')
        self.R=PolynomialRing(ZZ,[f'x_{i}_{j}' for i in range(k) for j in range(n)],order='lex')
        self.X=matrix(self.R,k,n,self.R.gens())
        self.minors={};self.polynomials={};self.reductions={}

    def polynomial(self,pair):
        if pair not in self.polynomials:
            for i in pair:
                if i not in self.minors:self.minors[i]=self.X.matrix_from_columns(self.columns[i]).det()
            poly=(self.minors[pair[0]]*self.minors[pair[1]]).dict()
            self.entries+=len(poly)
            if self.entries>self.limit:raise LimitExceeded('Quadratic straightening cache exceeds max_entries')
            self.polynomials[pair]=poly
        return self.polynomials[pair]

    def reduce(self,pair):
        if pair in self.index:return {self.index[pair]:ZZ.one()}
        if pair in self.reductions:return self.reductions[pair]
        work=dict(self.polynomial(pair));answer={};steps=0
        while work:
            exponent=max(work);coefficient=work[exponent]
            rows=[[j for j in range(self.n) for _ in range(exponent[i*self.n+j])] for i in range(self.k)]
            if any(len(row)!=2 for row in rows):raise ArithmeticError('Unexpected multidegree in straightening')
            a=self.column_index[tuple(row[0] for row in rows)]
            b=self.column_index[tuple(row[1] for row in rows)]
            mon=tuple(sorted((a,b)));index=self.index[mon]
            answer[index]=answer.get(index,0)+coefficient
            poly=self.polynomial(mon)
            if poly.get(exponent)!=1:raise ArithmeticError('Nonunit standard leading coefficient')
            for e,c in poly.items():
                value=work.get(e,0)-coefficient*c
                if value:work[e]=value
                else:work.pop(e,None)
            if work and max(work)>=exponent:raise ArithmeticError('Straightening did not decrease its leading term')
            steps+=1
            if steps>self.limit:raise LimitExceeded('Quadratic straightening step bound exceeded')
        answer={i:c for i,c in answer.items() if c}
        self.reductions[pair]=answer
        return answer

    def form_vector(self,B):
        result=vector(B.base_ring(),len(self.index))
        for i in range(B.nrows()):
            for j in range(i,B.ncols()):
                c=B[i,j] if i==j else B[i,j]+B[j,i]
                if c:
                    for index,a in self.reduce((i,j)).items():result[index]+=c*a
        return normalize(result)


def group_word_tree(rep,max_entries):
    """Enumerate small n-dimensional matrices; retain only words for averaging."""
    K=rep['field'];generators=rep['generators'];n=generators[0].nrows();order=rep['group_order']
    if order*n*n>max_entries:raise LimitExceeded('Group enumeration exceeds max_entries')
    one=identity_matrix(K,n);one.set_immutable()
    seen={one:0};elements=[one];tree=[None];characters=[K.one()]
    for parent,M in enumerate(elements):
        for j,S in enumerate(generators):
            T=S*M;T.set_immutable();value=rep['generator_characters'][j]*characters[parent]
            if T in seen:
                if characters[seen[T]]!=value:raise ArithmeticError('Scalar character does not descend to the matrix group')
                continue
            if len(elements)>=order:raise ArithmeticError('Generated matrix group exceeds its certified order')
            seen[T]=len(elements);elements.append(T);tree.append((parent,j));characters.append(value)
    if len(tree)!=order:raise ArithmeticError('Matrix group order differs from the character-search group')
    return tree,[1/c for c in characters]


def rational_components(V):
    """Separate a cyclotomic matrix into rational matrices (FLINT arithmetic).

    Combine proportional coefficient matrices; e.g. Q(sqrt(-7)) represented
    inside Q(zeta_7) often needs only two components rather than six.
    """
    K=V.base_ring();degree=K.degree();total=V.nrows()*V.ncols()
    arrays=[[QQ.zero()]*total for _ in range(degree)]
    for i,c in enumerate(V.list()):
        for j,a in enumerate(c.list()):arrays[j][i]=a
    grouped={}
    for j,values in enumerate(arrays):
        scale=next((a for a in values if a),None)
        if scale is None:continue
        M=matrix(QQ,V.nrows(),V.ncols(),values)/scale;M.set_immutable()
        grouped[M]=grouped.get(M,K.zero())+scale*K.gen()**j
    return [(c,M) for M,c in grouped.items() if c]


def reynolds_gram(V,weights):
    """Sum psi(g)^(-1) (g v)(g v)^t, without complex conjugation."""
    K=V.base_ring();trivial=all(w==1 for w in weights)
    if K is QQ:
        weighted=V if trivial else matrix(K,[w*row for w,row in zip(weights,V.rows())])
        return V.transpose()*weighted
    components=rational_components(V)
    B=matrix(K,V.ncols())
    if trivial:
        for i,(a,A) in enumerate(components):
            for j in range(i,len(components)):
                b,C=components[j];T=A.transpose()*C
                if i!=j:T=T+T.transpose()
                B+=(a*b)*T.change_ring(K)
    else:
        weighted=matrix(K,[w*row for w,row in zip(weights,V.rows())])
        weighted_components=rational_components(weighted)
        for a,A in components:
            for b,C in weighted_components:B+=(a*b)*(A.transpose()*C).change_ring(K)
    return B


class RationalOrbitAction:
    """Restrict the actual coefficient subfield to Q for fast orbit propagation.

    For example, a Q(sqrt(-7)) representation stored in Q(zeta_7) uses
    two rational blocks, not six. Closure under multiplication is checked
    exactly; no numerical field recognition is used.
    """
    def __init__(self,K,generators,max_entries):
        self.K=K;self.N=generators[0].nrows()
        if K is QQ:
            self.basis=[QQ.one()];self.generators=generators;return
        degree=K.degree();ambient=VectorSpace(QQ,degree)
        def coordinates(c):return vector(QQ,list(c.list())+[0]*(degree-len(c.list())))
        entries={K.one()}
        for A in generators:entries.update(A.list())
        W=ambient.span([coordinates(c) for c in entries])
        while True:
            elements=[sum(c*K.gen()**i for i,c in enumerate(row)) for row in W.basis()]
            extended=W+ambient.span([coordinates(a*b) for a in elements for b in elements])
            if extended.dimension()==W.dimension():break
            W=extended
        self.basis=elements;r=len(elements);N=self.N
        if elements[0]!=1:raise ArithmeticError('Coefficient-field basis must start with 1')
        if len(generators)*(r*N)**2>max_entries:raise LimitExceeded('Rational restriction matrices exceed max_entries')
        coordinate_cache={c:list(W.coordinate_vector(coordinates(c))) for c in entries}
        products=[[list(W.coordinate_vector(coordinates(a*b))) for b in elements] for a in elements]
        self.generators=[]
        for A in generators:
            arrays=[[QQ.zero()]*(N*N) for _ in range(r)]
            for pos,c in enumerate(A.list()):
                for j,a in enumerate(coordinate_cache[c]):arrays[j][pos]=a
            components=[matrix(QQ,N,N,values) for values in arrays]
            blocks=[]
            for i in range(r):
                row=[]
                for j in range(r):
                    block=matrix(QQ,N,N)
                    for t,C in enumerate(components):
                        if products[t][j][i]:block+=products[t][j][i]*C
                    row.append(block)
                blocks.append(row)
            self.generators.append(block_matrix(blocks))

    def components(self,v,tree):
        N=self.N;r=len(self.basis)
        if not hasattr(self,'integer_generators'):
            self.integer_generators=[]
            for A in self.generators:
                denominator=A.denominator()
                self.integer_generators.append(((denominator*A).change_ring(ZZ),denominator))
        initial_denominator=lcm([QQ(x).denominator() for x in v])
        initial=vector(ZZ,[initial_denominator*x for x in v]+[0]*(N*(r-1)))
        orbit=[initial];denominators=[initial_denominator]
        for parent,j in tree[1:]:
            A,denominator=self.integer_generators[j]
            orbit.append(A*orbit[parent]);denominators.append(denominator*denominators[parent])
        denominator=lcm(denominators)
        V=matrix(ZZ,[(denominator//d)*row for d,row in zip(denominators,orbit)])
        return [(a/denominator,V.matrix_from_columns(range(j*N,(j+1)*N))) for j,a in enumerate(self.basis)]

    def gram(self,v,tree,weights):
        components=self.components(v,tree);K=self.K
        if K is QQ:
            a,V=components[0]
            if all(w==1 for w in weights):B=V.transpose()*V
            else:B=V.transpose()*matrix(ZZ,[ZZ(w)*row for w,row in zip(weights,V.rows())])
            return (a*a)*B.change_ring(QQ)
        if all(w==1 for w in weights):
            B=matrix(K,self.N)
            for i,(a,A) in enumerate(components):
                for j in range(i,len(components)):
                    b,C=components[j];T=A.transpose()*C
                    if i!=j:T=T+T.transpose()
                    B+=(a*b)*T.change_ring(K)
            return B
        V=sum((a*C.change_ring(K) for a,C in components),matrix(K,len(tree),self.N))
        return reynolds_gram(V,weights)


def verify_form(B,generators,characters):
    K=B.base_ring()
    if B!=B.transpose():return False
    if K is QQ:return all(A*B*A.transpose()==chi*B for A,chi in zip(generators,characters))
    # Avoid generic number-field matrix multiplication in the certificate check.
    Bparts=rational_components(B)
    for A,chi in zip(generators,characters):
        Aparts=rational_components(A);image=matrix(K,B.nrows())
        for a,P in Aparts:
            for b,Q in Bparts:
                PQ=P*Q
                for c,T in Aparts:image+=(a*b*c)*(PQ*T.transpose()).change_ring(K)
        if image!=chi*B:return False
    return True

def invariant_quadrics(rep,k,n,multiplicity,*,max_sections=None,max_entries=2000000,seed=0,max_trials=0,verbose=False,return_certificate=False):
    start=time.perf_counter();K=rep['field'];dim=section_dimension(k,n,2)
    if max_entries<1 or max_trials<0:raise ValueError('max_entries must be positive and quadratic_trials nonnegative')
    if max_sections is not None and dim>max_sections:
        raise LimitExceeded(f'Section space has dimension {dim}; explicit max_sections is {max_sections}')
    if multiplicity<1:raise ValueError('A positive invariant multiplicity is required')
    if dim*multiplicity>max_entries:raise LimitExceeded('Stored quadratic invariant basis exceeds max_entries')
    columns,monomials=standard_monomials(k,n,2,dim);N=len(columns)
    degree=1 if K is QQ else K.degree()
    if max(rep['group_order']*N*degree,len(rep['generators'])*N*N,multiplicity*N*N)>max_entries:
        raise LimitExceeded('Reynolds orbit or exterior generator matrices exceed max_entries')
    tree,weights=group_word_tree(rep,max_entries)
    generators=[compound(M,k,columns) for M in rep['generators']]
    action=RationalOrbitAction(K,generators,max_entries)
    straightener=QuadraticStraightener(k,n,columns,monomials,max_entries)
    rng=random.Random(seed);points=[];plucker_points=[]
    for _ in range(max(8,multiplicity+4)):
        point=matrix(QQ,k,n,[rng.randrange(-2,3) for _ in range(k*n)])
        points.append(point);plucker_points.append(vector(K,minors(point,columns)))
    basis=[];evaluations=[];lifts=[];bound=max_trials or max(32,4*multiplicity+8)
    for trial in range(bound):
        v=vector(QQ,[rng.randrange(-2,3) for _ in range(N)])
        B=action.gram(v,tree,weights)
        values=vector(K,[(p*B).dot_product(p) for p in plucker_points])
        if matrix(K,evaluations+[values]).rank()==len(evaluations):continue
        # Exact ambient semi-invariance before quotienting by the Plucker ideal.
        if not verify_form(B,generators,rep['generator_characters']):
            raise ArithmeticError('Reynolds quadratic failed exact generator verification')
        reduced=straightener.form_vector(B)
        # Also certify that the integer straightening identities preserve all
        # evaluations used to prove quotient independence (up to normalization).
        raw=vector(K,[sum(c*p[i]*p[j] for (i,j),c in zip(monomials,reduced) if c) for p in plucker_points])
        nonzero=next(i for i,c in enumerate(values) if c);ratio=raw[nonzero]/values[nonzero]
        if not ratio or raw!=ratio*values:raise ArithmeticError('Quadratic quotient evaluation mismatch')
        basis.append(reduced);evaluations.append(values)
        if return_certificate:lifts.append(ratio*B)
        if verbose:print(f'Quadratic invariants: {len(basis)}/{multiplicity} certified',flush=True)
        if len(basis)==multiplicity:break
    if len(basis)!=multiplicity:
        raise LimitExceeded(f'Found {len(basis)}/{multiplicity} independent quadratic sections in {bound} trials; increase quadratic_trials or change seed')
    # RREF is canonical for the complete quotient invariant space, independent
    # of random lifts. It also tends to remove the dense Reynolds coefficients.
    original_basis=matrix(K,basis)
    basis=[normalize(v) for v in original_basis.echelon_form().rows() if v]
    evidence={'method':'exact character-weighted Reynolds forms and Plucker straightening',
              'group_order':len(tree),'plucker_dimension':N,'ambient_quadratic_dimension':N*(N+1)//2,
              'section_dimension':dim,'invariant_dimension':len(basis),'trials':trial+1,'seed':seed,
              'coefficient_subfield_degree':len(action.basis),'evaluation_rank':len(evaluations),'evaluation_points':[encode_matrix(p) for p in points],
              'verification':'Exact ambient generator identities; integer generic-minor straightening; exact rank on decomposable evaluation points equals character multiplicity',
              'straightening_relations':len(straightener.reductions),'seconds':time.perf_counter()-start}
    result=(columns,monomials,basis,evidence)
    if return_certificate:
        return result+({'original_basis':original_basis,'lifts':lifts,'straightener':straightener,'generators':generators},)
    return result


def construct_quadratic(candidate,args):
    from .representations import realize
    from .models import make_model
    from .construct import finish_model
    k,n=candidate['k'],candidate['n'];verbose=getattr(args,'verbose',False)
    dim=section_dimension(k,n,2)
    if args.max_sections is not None and dim>args.max_sections:
        raise LimitExceeded(f'Section space has dimension {dim}; explicit max_sections is {args.max_sections}')
    if dim*candidate['multiplicity']>args.max_entries:
        raise LimitExceeded('Stored quadratic invariant basis exceeds max_entries')
    print('Constructing quadratic matrices for',candidate['group'],flush=True)
    rep=realize(candidate,args.gap_group);K=rep['field']
    columns,monomials,basis,evidence,certificate=invariant_quadrics(rep,k,n,candidate['multiplicity'],
        max_sections=args.max_sections,max_entries=args.max_entries,seed=args.seed,
        max_trials=getattr(args,'quadratic_trials',0),verbose=verbose,return_certificate=True)
    v=choose_vector(basis,args.seed,args.trials,args.coefficients.split(',') if args.coefficients else None,args.random_section)
    # Recover an exactly verified ambient lift of the selected quotient vector.
    # This avoids re-expanding a large transformed polynomial on a dense chart.
    E=certificate['original_basis'];pivots=E.pivots()
    coefficients=E.matrix_from_columns(pivots).transpose().solve_right(vector(K,[v[j] for j in pivots]))
    if coefficients*E!=v:raise ArithmeticError('Selected quadric is outside the certified invariant space')
    lift=sum((c*B for c,B in zip(coefficients,certificate['lifts'])),matrix(K,len(columns)))
    if certificate['straightener'].form_vector(lift)!=normalize(v):raise ArithmeticError('Selected quadratic lift failed straightening')
    if not verify_form(lift,certificate['generators'],rep['generator_characters']):raise ArithmeticError('Selected quadratic lift failed generator verification')
    original_v=v
    C=identity_matrix(K,n);matrices=rep['generators'];label='quadratic-Reynolds + standard-monomial RREF'
    if not args.coefficients and not args.random_section:
        simplified=diagonal_simplification(K,columns,monomials,v,basis)
        if simplified is not None:
            C,v,basis=simplified;Ci=C.inverse();matrices=[Ci*M*C for M in matrices];label+=' + diagonal rescaling'
    recovered=vector(K,[c*prod(C[j,j] for t in mon for j in columns[t]) for c,mon in zip(v,monomials)])
    if normalize(recovered)!=normalize(original_v):raise ArithmeticError('Diagonal coordinate-change certificate failed')
    evidence['selected_equation_verification']='Exact semi-invariant ambient lift; integer Plucker straightening; exact diagonal coordinate-change identity'
    model=make_model(candidate,rep,label,K,C,columns,monomials,v,basis,[])
    model['quadratic_construction']=evidence
    return finish_model(model,args)
