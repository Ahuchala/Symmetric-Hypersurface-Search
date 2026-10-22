"""Conservative chart reduction using exactly verified monomial symmetries."""
import itertools
from sage.all import Graph
from .common import decode_matrix, encode
from .models import coefficients
from .chart_engine import coefficient_reduction


def sign_sort(indices):
    return (-1)**sum(indices[i]>indices[j] for i in range(len(indices)) for j in range(i+1,len(indices)))


def verified_action(model,perm,weights,prime):
    """Sufficient ambient-polynomial identity; no equality modulo relations assumed."""
    n=model['n']
    if sorted(perm)!=list(range(n)) or len(weights)!=n or not all(weights):return None
    perm=[int(i) for i in perm]
    K,terms=coefficients(model)
    # Reduce both the equation and the proposed diagonal weights at one embedding.
    probe=dict(model);probe['terms']=list(model['terms'])+[
        {'monomial':[0],'coefficient':encode(w)} for w in weights]
    F,convert,embedding=coefficient_reduction(probe,prime)
    weights=[convert(w) for w in weights]
    if not all(weights):return None
    cols=model['plucker_indices'];lookup={tuple(I):j for j,I in enumerate(cols)}
    mapping=[]
    for I in cols:
        image=[perm[i] for i in I];scale=F(sign_sort(image))
        for i in I:scale*=weights[i]
        mapping.append((lookup[tuple(sorted(image))],scale))
    original={};transformed={}
    for mon,c in terms:
        c=convert(c);key=tuple(sorted(mon))
        original[key]=original.get(key,F.zero())+c
        target=[]
        for j in mon:
            index,scale=mapping[j];target.append(index);c*=scale
        key=tuple(sorted(target));transformed[key]=transformed.get(key,F.zero())+c
    original={k:v for k,v in original.items() if v};transformed={k:v for k,v in transformed.items() if v}
    if not original or original.keys()!=transformed.keys():return None
    first=next(iter(original));scalar=transformed[first]/original[first]
    if not scalar or any(transformed[k]!=scalar*c for k,c in original.items()):return None
    return {'permutation':list(perm),'weights':[str(w) for w in weights],
            'scalar':str(scalar),'embedding':embedding,
            'verification':'exact ambient polynomial semi-invariance over the checking field'}


def chart_orbits(model,prime,enabled=False,permutations=None):
    pivots=list(itertools.combinations(range(model['n']),model['k']))
    if not enabled and not permutations:return [[I] for I in pivots],[]
    K,terms=coefficients(model);n=model['n'];candidates=[]
    for perm in permutations or []:
        proof=verified_action(model,list(perm),[K.one()]*n,prime)
        if proof is None:raise ValueError('Supplied coordinate permutation is not a verified symmetry')
        candidates.append(proof)
    if enabled:
        # Saved generators need not be monomial. Ignore those that do not
        # permute standard charts, rather than assuming all group elements do.
        for encoded in model.get('generators',[]):
            M=decode_matrix(encoded,K);perm=[];weights=[]
            for j in range(n):
                entries=[i for i in range(n) if M[i,j]]
                if len(entries)!=1:break
                perm.append(entries[0]);weights.append(M[entries[0],j])
            if len(perm)==n:
                try:proof=verified_action(model,perm,weights,prime)
                except (ValueError,TypeError,ArithmeticError):proof=None
                if proof:candidates.append(proof)
        # Bounded-size support-graph discovery for sparse linear sections.
        # Graph automorphisms are proposals only: every one is checked above.
        if model['d']==1 and n<=12 and len(terms)<=100:
            edges=[model['plucker_indices'][mon[0]] for mon,c in terms if c]
            if edges:
                G=Graph();G.add_vertices(range(n+len(edges)))
                G.add_edges((v,n+j) for j,I in enumerate(edges) for v in I)
                group=G.automorphism_group(partition=[list(range(n)),list(range(n,n+len(edges)))])
                for g in group.gens():
                    perm=[int(g(i)) for i in range(n)]
                    proof=verified_action(model,perm,[K.one()]*n,prime)
                    if proof:candidates.append(proof)
    unique={tuple(c['permutation']):c for c in candidates if c['permutation']!=list(range(n))}
    proofs=list(unique.values());unseen=set(pivots);orbits=[]
    while unseen:
        seed=min(unseen);orbit={seed};todo=[seed]
        while todo:
            I=todo.pop()
            for proof in proofs:
                J=tuple(sorted(proof['permutation'][i] for i in I))
                if J not in orbit:orbit.add(J);todo.append(J)
        unseen-=orbit;orbits.append(sorted(orbit))
    return orbits,proofs
