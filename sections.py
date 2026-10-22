"""Standard monomials and exact semi-invariants in the Grassmannian ring."""
import itertools
import random
from sage.all import QQ, matrix, vector, PolynomialRing, identity_matrix, prod, binomial, gcd, lcm
from .common import LimitExceeded, conductor, field_from


def section_dimension(k,n,d):
    return int(prod(QQ(n+j-i)/QQ(d-j+k-i+1) for i in range(1,k+1) for j in range(1,d+1)))


def standard_monomials(k,n,d,max_sections=1000):
    dim=section_dimension(k,n,d)
    if dim>max_sections:
        raise LimitExceeded(f'Section space has dimension {dim}; construction limit is {max_sections}. Character searches remain available; raise --max-sections to construct.')
    columns=list(itertools.combinations(range(n),k))
    if d==1:
        return columns,[(i,) for i in range(len(columns))]
    out=[]
    def visit(prefix):
        if len(prefix)==d:
            out.append(tuple(prefix));return
        start=prefix[-1] if prefix else 0
        for j in range(start,len(columns)):
            if not prefix or all(a<=b for a,b in zip(columns[prefix[-1]],columns[j])):
                visit(prefix+[j])
    visit([])
    assert len(out)==dim
    return columns,out


def chart_matrix(K,k,n,pivot):
    pivot=tuple(pivot)
    free=[j for j in range(n) if j not in pivot]
    R=PolynomialRing(K,[f'x_{i}_{j}' for i in range(k) for j in free],order='degrevlex')
    X=matrix(R,k,n)
    for i,j in enumerate(pivot):X[i,j]=1
    for x,(i,j) in zip(R.gens(),itertools.product(range(k),free)):X[i,j]=x
    return R,X


def minors(X,columns):
    rows=list(range(X.nrows()))
    return [X.matrix_from_rows_and_columns(rows,J).det() for J in columns]


def compound(M,k,columns):
    return matrix(M.base_ring(),[[M.matrix_from_rows_and_columns(I,J).det() for J in columns] for I in columns])


def normalize(v):
    first=next((x for x in v if x),None)
    if first is None:return v
    v=v/first
    if all(x in QQ for x in v):
        denom=lcm([QQ(x).denominator() for x in v]);v=denom*v
        content=gcd([QQ(x).numerator() for x in v if x]);v=v/content
    return v


def invariant_space(K,matrices,characters,k,n,d,columns,monomials,allowed=None,max_entries=2000000):
    allowed=list(range(len(monomials))) if allowed is None else allowed
    if not allowed:
        return []
    if d==1:
        blocks=[]
        for M,p in zip(matrices,characters):
            A=compound(M,k,columns)-p
            blocks.append(A.matrix_from_columns(allowed))
        equations=blocks[0]
        for b in blocks[1:]:equations=equations.stack(b)
    else:
        R,X=chart_matrix(K,k,n,range(k));p0=minors(X,columns)
        basis=[prod(p0[j] for j in monomials[i]) for i in allowed]
        entries={};offset=0
        for M,char in zip(matrices,characters):
            images=minors(X*M,columns)
            polys=[prod(images[j] for j in monomials[i])-char*b for i,b in zip(allowed,basis)]
            exponents=sorted(set().union(*(f.dict() for f in polys)))
            if (offset+len(exponents))*len(allowed)>max_entries:
                raise LimitExceeded('Expanded invariance system exceeds --max-entries')
            lookup={e:offset+j for j,e in enumerate(exponents)}
            for col,f in enumerate(polys):
                entries.update({(lookup[e],col):c for e,c in f.dict().items()})
            offset+=len(exponents)
        equations=matrix(K,offset,len(allowed),entries,sparse=True)
    result=[]
    for v in equations.right_kernel().basis():
        full=vector(K,len(monomials))
        for j,x in zip(allowed,v):full[j]=x
        result.append(normalize(full))
    return result


def normalizer_basis(rep, cl, L, attempts=250):
    """Bounded exact search for a cyclic normalizer; preserve eigenspace orbits."""
    g=cl['matrix'];order=cl['order'];K=g.base_ring();n=g.nrows()
    powers={a:g**a for a in range(2,order) if gcd(a,order)==1}
    if not powers:return None
    h=identity_matrix(K,n);rng=random.Random(0)
    found=None
    for _ in range(attempts):
        h=h*rep['generators'][rng.randrange(len(rep['generators']))]
        hg=h*g
        for a,ga in powers.items():
            if hg!=ga*h:continue
            length=1;power=a%order
            while power!=1:
                length+=1;power=power*a%order
            if h**length==1:
                found=(a,h.change_ring(L));break
        if found:break
    if found is None:return None
    a,h=found;inverse=int(pow(a,-1,order));root=L.zeta(order);g=g.change_ring(L)
    spaces={j:list((g-root**j).right_kernel().basis()) for j in range(order)}
    remaining={j for j in spaces if spaces[j]};vectors=[];weights=[]
    while remaining:
        j=min(remaining);orbit=[j];nextj=inverse*j%order
        while nextj!=j:
            orbit.append(nextj);nextj=inverse*nextj%order
        for v in spaces[j]:
            for t,weight in enumerate(orbit):
                vectors.append(h**t*v);weights.append(root**weight)
        remaining.difference_update(orbit)
    C=matrix(L,vectors).transpose()
    if not C.det():return None
    return C,weights


def candidate_bases(rep,k,n,d,columns,monomials,tries,max_conductor):
    K=rep['field']
    yield 'original',K,identity_matrix(K,n),None
    for label,C in rep['seed_bases']:
        yield label,K,C,None
    if tries <= 0:
        return
    ordered=sorted((c for c in rep['classes'] if c['order']>1),key=lambda c:-c['order'])
    seen=set();used=0
    for cl in ordered:
        order=cl['order']
        if order in seen:continue
        seen.add(order)
        m=int(lcm(conductor(K),order))
        if m>max_conductor:continue
        L=field_from(m);g=cl['matrix'].change_ring(L)
        root=L.zeta(order)
        vectors=[];weights=[]
        for j in range(order):
            eig=(g-root**j).right_kernel().basis()
            vectors.extend(eig);weights.extend([root**j]*len(eig))
        if len(vectors)!=n:continue
        C=matrix(L,vectors).transpose()
        pw=[prod(weights[j] for j in col) for col in columns]
        wanted=L(cl['character'])
        allowed=[i for i,mon in enumerate(monomials) if prod(pw[j] for j in mon)==wanted]
        yield f'order-{order}-eigenbasis',L,C,allowed
        orbit_basis=normalizer_basis(rep,cl,L)
        if orbit_basis is not None:
            C,weights=orbit_basis
            pw=[prod(weights[j] for j in col) for col in columns]
            allowed=[i for i,mon in enumerate(monomials) if prod(pw[j] for j in mon)==wanted]
            yield f'order-{order}-normalizer-basis',L,C,allowed
        used+=1
        if used>=tries:return


def choose_vector(basis,seed,trials,coefficients=None,random_section=False):
    K=basis[0].base_ring();rng=random.Random(seed)
    if coefficients is not None:
        if len(coefficients)!=len(basis):raise ValueError(f'Need {len(basis)} coefficients')
        out=sum((K(c)*v for c,v in zip(coefficients,basis)),vector(K,len(basis[0])))
        if not out:raise ValueError('The requested combination is zero')
        return normalize(out)
    if random_section:
        coeffs=[rng.randint(-3,3) for _ in basis]
        if not any(coeffs):coeffs[0]=1
        return normalize(sum((c*v for c,v in zip(coeffs,basis)),vector(K,len(basis[0]))))
    choices=list(basis)
    for _ in range(trials if len(basis)>1 else 0):
        coeffs=[rng.randint(-2,2) for _ in basis]
        v=sum((c*v for c,v in zip(coeffs,basis)),vector(K,len(basis[0])))
        if v:choices.append(normalize(v))
    return min(choices,key=score)


def score(v):
    def height(x):
        coefficients=[QQ(x)] if x in QQ else x.list()
        return sum(abs(int(c.numerator())).bit_length()+int(c.denominator()).bit_length() for c in coefficients if c)+(0 if x in QQ else 4*len(coefficients))
    return (sum(bool(x) for x in v),sum(height(x) for x in v if x))


def verify_section(K,matrices,characters,k,n,d,columns,monomials,v):
    if d==1:
        for M,char in zip(matrices,characters):
            assert compound(M,k,columns)*v==char*v
    else:
        _,X=chart_matrix(K,k,n,range(k));p0=minors(X,columns)
        f=sum(c*prod(p0[j] for j in mon) for mon,c in zip(monomials,v) if c)
        assert f
        for M,char in zip(matrices,characters):
            p=minors(X*M,columns)
            assert sum(c*prod(p[j] for j in mon) for mon,c in zip(monomials,v) if c)==char*f


def diagonal_simplification(K,columns,monomials,v,basis):
    """Exact multiplicative elimination using unit pivots; no numerical roots.

    Try to make nonzero coefficients equal up to sign by rescaling the original
    n coordinates. Unsolvable residual equations are retained: the caller accepts
    the change only if it improves the actual verified coefficient score.
    """
    from sage.all import diagonal_matrix, ZZ
    n=1+max(max(I) for I in columns)
    support=[i for i,c in enumerate(v) if c]
    def exponents(mon):
        out=[0]*n
        for j in mon:
            for t in columns[j]:out[t]+=1
        return out
    powers=[exponents(mon) for mon in monomials]
    def signed(c):
        vals=[QQ(c)] if c in QQ else c.list()
        return c if next(a for a in vals if a)>0 else -c
    first=support[0];base=signed(v[first])
    rows=[[a-b for a,b in zip(powers[i],powers[first])] for i in support[1:]]
    rhs=[signed(v[i])/base for i in support[1:]]
    row=0;pivots=[]
    while row<len(rows):
        choice=next(((i,j) for i in range(row,len(rows)) for j in range(n) if abs(rows[i][j])==1 and j not in pivots),None)
        if choice is None:break
        i,j=choice;rows[row],rows[i]=rows[i],rows[row];rhs[row],rhs[i]=rhs[i],rhs[row]
        if rows[row][j]==-1:
            rows[row]=[-a for a in rows[row]];rhs[row]=1/rhs[row]
        for i in range(len(rows)):
            if i==row:continue
            a=rows[i][j]
            if a:
                rows[i]=[x-a*y for x,y in zip(rows[i],rows[row])]
                rhs[i]/=rhs[row]**a
        pivots.append(j);row+=1
    scales=[K(1)]*n
    for i,j in enumerate(pivots):scales[j]=rhs[i]
    divisors=[prod(scales[j]**a for j,a in enumerate(exp)) for exp in powers]
    new=normalize(vector(K,[c/t for c,t in zip(v,divisors)]))
    if score(new)>=score(v):return None
    newbasis=[normalize(vector(K,[c/t for c,t in zip(u,divisors)])) for u in basis]
    return diagonal_matrix(K,scales),new,newbasis
