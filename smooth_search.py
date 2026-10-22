"""Bounded discovery of certified smooth hypersurfaces with faithful action."""
import inspect
import random
from .common import decode, encode


def combination(X, coefficients):
    """Use the certified basis in its saved coordinates; discard old checks."""
    from .api import Hypersurface
    model=X.data
    if len(coefficients)!=len(model['invariant_basis']):
        raise ValueError('One coefficient per invariant basis element is required')
    terms={}
    for a,basis in zip(coefficients,model['invariant_basis']):
        for term in basis:
            mon=tuple(term['monomial'])
            terms[mon]=terms.get(mon,X.field.zero())+X.field(a)*decode(term['coefficient'],X.field)
    terms={m:c for m,c in terms.items() if c}
    if not terms:raise ValueError('Zero section')
    scale=terms[min(terms)]
    model['terms']=[{'monomial':list(m),'coefficient':encode(c/scale)} for m,c in sorted(terms.items())]
    model.pop('faithfulness_evidence',None)
    model['verification']['smoothness']='unchecked'
    model['verification']['faithfulness_on_hypersurface']='unchecked'
    model['section_selection']={'method':'combination of certified saved invariant basis',
                                'coefficients':[encode(X.field(c)) for c in coefficients]}
    if 'quadratic_construction' in model:
        model['quadratic_construction']['selected_equation_verification']='Linear combination of certified invariant basis in saved coordinates'
    return Hypersurface(model)


def equation_key(X):
    terms=sorted((tuple(t['monomial']),decode(t['coefficient'],X.field)) for t in X._model['terms'])
    scale=terms[0][1]
    return tuple((m,tuple(encode(c/scale))) for m,c in terms)


def find_smooth(results, *, sections_per_candidate, max_candidates, seed, prime,
                seconds, max_charts, algorithm, symmetry, workers,
                faithfulness_charts, construction_options, verbose):
    from .api import Candidate
    sections_per_candidate=int(sections_per_candidate)
    max_candidates=len(results) if max_candidates is None else int(max_candidates)
    if sections_per_candidate<1 or max_candidates<0:
        raise ValueError('sections_per_candidate must be positive; max_candidates nonnegative')
    if seconds<=0 or max_charts<0 or faithfulness_charts<0 or workers<1:
        raise ValueError('Invalid checking budget')
    if algorithm not in ('auto','std','slimgb'):raise ValueError('Unknown algorithm')
    options=dict(construction_options or {})
    reserved={'seed','coefficients','random_section','faithfulness_charts','verbose'}
    if reserved.intersection(options):raise ValueError('Construction options reserved by find_smooth: '+', '.join(sorted(reserved.intersection(options))))
    inspect.signature(Candidate.construct).bind(None,**options)
    options.update(seed=int(seed),faithfulness_charts=0,verbose=False)
    report={'status':'running','search_complete':results.complete,'seed':int(seed),
            'sections_per_candidate':sections_per_candidate,'max_candidates':max_candidates,
            'checking':{'prime':int(prime),'seconds_per_chart':float(seconds),'max_charts':int(max_charts),
                        'algorithm':algorithm,'symmetry':bool(symmetry),'workers':int(workers),
                        'faithfulness_charts':int(faithfulness_charts)},'attempts':[]}
    results.last_smooth_search=report
    rng=random.Random(int(seed))
    try:
        for index,candidate in enumerate(results[:max_candidates]):
            if verbose:print(f'[{index}] {candidate.group}: constructing invariant space',flush=True)
            try:base=candidate.construct(**options)
            except (RuntimeError,ValueError) as exc:
                report['attempts'].append({'candidate_index':index,'candidate_id':candidate.id,
                    'group':candidate.group,'stage':'construction','status':'inconclusive','error':str(exc)})
                continue
            seen=set();dimension=len(base._model['invariant_basis'])
            for trial in range(1 if dimension==1 else sections_per_candidate):
                coeffs=None
                if trial==0:X=base
                else:
                    for _ in range(100):
                        coeffs=[rng.randint(-5,5) for _ in range(dimension)]
                        if not any(coeffs):continue
                        X=combination(base,coeffs)
                        if equation_key(X) not in seen:break
                    else:break
                seen.add(equation_key(X))
                attempt={'candidate_index':index,'candidate_id':candidate.id,'group':candidate.group,
                         'section':trial,'coefficients':coeffs,'status':'checking'}
                report['attempts'].append(attempt)
                try:
                    faithful=X.check_faithfulness(max_charts=faithfulness_charts)
                    attempt['faithfulness']=faithful.report
                    if faithful.is_faithful is not True:
                        attempt['status']='not_faithful' if faithful.is_faithful is False else 'inconclusive'
                    else:
                        smooth=X.check_smoothness(prime=prime,seconds=seconds,max_charts=max_charts,
                            algorithm=algorithm,symmetry=symmetry,workers=workers)
                        attempt['smoothness']=smooth.report
                        attempt['status']=smooth.status
                        if smooth.is_smooth is True:
                            attempt['status']='smooth_and_faithful';report['status']='found'
                            if verbose:print(f'[{index}] section {trial}: smooth and faithful',flush=True)
                            return X
                except (RuntimeError,ValueError) as exc:
                    attempt['status']='inconclusive';attempt['error']=str(exc)
                if verbose:print(f'[{index}] section {trial}: {attempt["status"]}',flush=True)
        report['status']='no_certificate_within_bounds'
        return None
    except BaseException:
        report['status']='interrupted'
        raise
