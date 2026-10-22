import json
from sage.all import vector
from .common import write_json
from .representations import realize
from .sections import standard_monomials, candidate_bases, invariant_space, choose_vector, score, verify_section, diagonal_simplification
from .models import make_model,print_equation
from .faithfulness import check_faithfulness


def construct(args):
    report=json.loads(open(args.search).read())
    matches=[c for c in report['candidates'] if c['id']==args.candidate]
    if len(matches)!=1:raise ValueError('Select exactly one candidate ID from the search file')
    return construct_candidate(matches[0],args)


def construct_candidate(candidate,args):
    k,n,d=(candidate[x] for x in ('k','n','d'))
    method=getattr(args,'quadratic_method','auto')
    if method not in ('auto','reynolds','expanded'):raise ValueError('Unknown quadratic_method')
    if method=='reynolds' and d!=2:raise ValueError('Reynolds quadratic construction requires d=2')
    if d==2 and method!='expanded':
        from .quadratics import construct_quadratic
        return construct_quadratic(candidate,args)
    columns,monomials=standard_monomials(k,n,d,1000 if args.max_sections is None else args.max_sections)
    print('Constructing matrices for',candidate['group'],flush=True)
    rep=realize(candidate,args.gap_group)
    attempts=[];best=None
    for label,K,C,allowed in candidate_bases(rep,k,n,d,columns,monomials,args.basis_tries,args.max_conductor):
        if getattr(args,'verbose',False):
            print(f'Trying {label}; {len(monomials) if allowed is None else len(allowed)} candidate monomials',flush=True)
        try:
            Ci=C.inverse();matrices=[Ci*M.change_ring(K)*C for M in rep['generators']]
            chars=[K(c) for c in rep['generator_characters']]
            basis=invariant_space(K,matrices,chars,k,n,d,columns,monomials,allowed,args.max_entries)
            if len(basis)!=candidate['multiplicity']:
                raise ArithmeticError(f'Character screen predicted {candidate["multiplicity"]} sections, explicit calculation found {len(basis)}')
            coeffs=args.coefficients.split(',') if args.coefficients else None
            v=choose_vector(basis,args.seed,args.trials,coeffs,args.random_section)
            simplified=None if args.coefficients or args.random_section else diagonal_simplification(K,columns,monomials,v,basis)
            if simplified is not None:
                S,v,basis=simplified
                C=C*S;Si=S.inverse()
                matrices=[Si*M*S for M in matrices]
                label+=' + diagonal rescaling'
            verify_section(K,matrices,chars,k,n,d,columns,monomials,v)
            current_score=score(v)
            attempts.append({'basis':label,'status':'verified','nonzero_terms':current_score[0],'coefficient_complexity':current_score[1]})
            if best is None or current_score<best[0]:
                best=(current_score,label,K,C,v,basis)
            # Explicit combinations are tied to one basis; never silently change them.
            if args.coefficients or args.random_section:break
        except (RuntimeError,ValueError) as exc:
            attempts.append({'basis':label,'status':'not_completed','error':str(exc)})
            print(' ',exc,flush=True)
    if best is None:raise RuntimeError('No construction completed: '+str(attempts))
    _,label,K,C,v,basis=best
    model=make_model(candidate,rep,label,K,C,columns,monomials,v,basis,attempts)
    return finish_model(model,args)


def finish_model(model,args):
    if args.output:
        write_json(args.output,model)
    if args.faithfulness_charts:
        print('Checking faithfulness on the reduced hypersurface',flush=True)
        evidence=check_faithfulness(model,args.faithfulness_charts)
        model['faithfulness_evidence']=evidence
        model['verification']['faithfulness_on_hypersurface']=evidence['status']
        if args.output:
            write_json(args.output,model)
    print_equation(model)
    if args.output:
        print('Saved equation, basis, matrices and verification to',args.output)
    return model
