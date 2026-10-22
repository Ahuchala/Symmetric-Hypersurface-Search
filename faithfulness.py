"""Pointwise kernel on the reduced hypersurface, not just the ambient space."""
import itertools
from sage.all import gcd, prod
from .common import field_from, decode_matrix
from .models import chart_polynomial
from .sections import minors


def check_faithfulness(model,max_charts=1):
    if not model.get('classes'):
        return {'status':'unchecked','reason':'No group conjugacy-class matrices supplied'}
    K=field_from(model['field_conductor']);n=model['n'];k=model['k']
    classes=[(i,cl,decode_matrix(cl['matrix'],K)) for i,cl in enumerate(model['classes']) if cl['order']!=1]
    pending={i:(cl,M) for i,cl,M in classes};witnesses=[];checked=0
    for pivot in itertools.combinations(range(n),k):
        if checked>=max_charts:break
        checked+=1
        R,X,p,f=chart_polynomial(model,pivot)
        if not f:raise ValueError('The proposed equation vanishes identically on the Grassmannian')
        if f.is_constant():continue
        repeated=f
        for x in R.gens():repeated=repeated.gcd(f.derivative(x))
        reduced=f//repeated
        base=model['plucker_indices'].index(list(pivot))
        for i,(cl,M) in list(pending.items()):
            images=minors(X*M,model['plucker_indices'])
            for j in range(len(p)):
                # p_base=1: fixed points satisfy g(p_j)=p_j*g(p_base).
                difference=images[j]-p[j]*images[base]
                remainder=difference % reduced if R.ngens()==1 else difference.reduce([reduced])
                if remainder:
                    witnesses.append({'class_index':i,'pivot':list(pivot),'plucker_index':j,'nonzero_remainder_terms':len(remainder.dict())})
                    del pending[i]
                    break
        if not pending:
            return {'status':'faithful','method':'Exact characteristic-zero non-fixed-locus witness for every nonidentity conjugacy class on the reduced hypersurface','charts_checked':checked,'witnesses':witnesses,'kernel_order':1}
    total=sum(1 for _ in itertools.combinations(range(n),k))
    if checked==total:
        return {'status':'not_faithful','method':'Every affine chart checked on the reduced hypersurface','kernel_order':1+sum(cl['size'] for cl,M in pending.values()),'kernel_class_indices':list(pending),'charts_checked':checked,'witnesses':witnesses}
    return {'status':'inconclusive','reason':'Some classes lack a non-fixed-locus witness; increase --max-charts','unresolved_classes':list(pending),'charts_checked':checked,'witnesses':witnesses}
