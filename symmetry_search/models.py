"""Serializable equations with explicit Plucker indices and matrix provenance."""
import json
from sage.all import prod, PolynomialRing, matrix, vector, QQ
from .common import field_from, encode, decode, encode_matrix, decode_matrix, conductor
from .sections import chart_matrix, minors


def make_model(candidate,rep,label,K,C,columns,monomials,v,basis,attempts):
    Ci=C.inverse()
    return {
        'schema':1,'candidate':candidate,'k':candidate['k'],'n':candidate['n'],'d':candidate['d'],
        'field_conductor':conductor(K),'coordinate_convention':'Generator columns act on the n-dimensional coordinate representation W. Plucker columns are increasing k-subsets of 0,...,n-1. The geometric action is the dual action.',
        'basis_strategy':label,'basis_attempts':attempts,
        'change_of_basis':encode_matrix(C),
        'plucker_indices':[list(I) for I in columns],
        'terms':[{'monomial':list(mon),'coefficient':encode(c)} for mon,c in zip(monomials,v) if c],
        'invariant_basis':[[{'monomial':list(mon),'coefficient':encode(c)} for mon,c in zip(monomials,u) if c] for u in basis],
        'generators':[encode_matrix(Ci*M.change_ring(K)*C) for M in rep['generators']],
        'generator_characters':[encode(K(c)) for c in rep['generator_characters']],
        'classes':[{'matrix':encode_matrix(Ci*cl['matrix'].change_ring(K)*C),'size':cl['size'],'order':cl['order']} for cl in rep['classes']],
        'verification':{'nonzero_section':True,'semi_invariance':'exact_characteristic_zero','invariant_dimension':len(basis),'smoothness':'unchecked','faithfulness_on_hypersurface':'unchecked'},
    }


def read_model(path):
    data=json.loads(open(path).read())
    if data.get('schema')!=1 or 'terms' not in data:
        raise ValueError('Expected a symmetry-search equation file')
    return data


def coefficients(model):
    K=field_from(model['field_conductor'])
    return K,[(tuple(t['monomial']),decode(t['coefficient'],K)) for t in model['terms']]


def chart_polynomial(model,pivot,K=None,map_coefficient=None,*,sparse=False,terms=None):
    if terms is None:
        L,terms=coefficients(model)
        K=L if K is None else K
    elif K is None:
        raise ValueError("Prepared terms require a coefficient field")
    convert=K if map_coefficient is None else map_coefficient
    R,X=chart_matrix(K,model['k'],model['n'],pivot)
    if sparse:
        used=sorted({j for mon,c in terms if c for j in mon})
        values=minors(X,[model['plucker_indices'][j] for j in used])
        p=[None]*len(model['plucker_indices'])
        for j,value in zip(used,values):p[j]=value
    else:
        p=minors(X,model['plucker_indices'])
    f=sum((convert(c)*prod(p[j] for j in mon) for mon,c in terms if c),R.zero())
    return R,X,p,f


def print_equation(model):
    K,terms=coefficients(model)
    names=['p_'+ '_'.join(map(str,I)) for I in model['plucker_indices']]
    P=PolynomialRing(K,names,order='degrevlex')
    f=sum(c*prod(P.gen(j) for j in mon) for mon,c in terms)
    print(f'Gr({model["k"]},{model["n"]}), degree {model["d"]}; {len(terms)} monomials; basis: {model.get("basis_strategy","supplied")}')
    if any(c not in QQ for _,c in terms):
        print(f'Coefficient field: Q(z), with z a primitive {model["field_conductor"]}-th root of unity')
    print('f =',f)
    print('Semi-invariance:',model.get('verification',{}).get('semi_invariance','unchecked'))
    print('Faithfulness on the reduced hypersurface:',model.get('verification',{}).get('faithfulness_on_hypersurface','unchecked'))
    return str(f)


def latex_equation(model):
    from collections import Counter
    from sage.all import latex, QQ
    K,terms=coefficients(model)
    chunks=[]
    for mon,c in terms:
        negative=c in QQ and QQ(c)<0
        a=-c if negative else c
        factors=[]
        for j,e in sorted(Counter(mon).items()):
            symbol='p_{'+','.join(map(str,model['plucker_indices'][j]))+'}'
            factors.append(symbol+('^{'+str(e)+'}' if e!=1 else ''))
        coefficient='' if a==1 else (latex(a) if a in QQ else r'\left('+latex(a)+r'\right)')
        chunks.append((' - ' if negative else (' + ' if chunks else ''))+coefficient+''.join(factors))
    return 'f = '+''.join(chunks).lstrip()
