"""Bounded character search, including semi-invariants and projective kernels."""
import hashlib
import itertools
import json
import re
from sage.all import libgap
from .common import ROOT, LimitExceeded, write_json
from .obstructions import smooth_character_obstruction

_INITIALIZED = False

def initialize():
    global _INITIALIZED
    if not _INITIALIZED:
        if not bool(libgap.LoadPackage('ctbllib')):
            raise RuntimeError('GAP CTblLib is required')
        libgap.Read(str(ROOT/'symmetry_search/characters.g'))
        _INITIALIZED = True


def validate(k,n,d):
    if not (0 < k < n and d > 0):
        raise ValueError('Require 0 < k < n and d >= 1')


def constituents(degrees, n, reducible, limit):
    count = 0
    if not reducible:
        for i,degree in enumerate(degrees):
            if degree == n:
                yield [i+1]
        return
    eligible = [(i+1,degree) for i,degree in enumerate(degrees) if degree <= n]
    def visit(start, remaining, selected):
        nonlocal count
        if not remaining:
            count += 1
            if count > limit:
                raise LimitExceeded(f'More than {limit} degree-{n} characters; raise --max-representations')
            yield selected
            return
        for j in range(start,len(eligible)):
            i,degree = eligible[j]
            if degree <= remaining:
                yield from visit(j,remaining-degree,selected+[i])
    yield from visit(0,n,[])


def candidates_for_table(tbl, group_spec, k,n,d, reducible=False, max_representations=10000,
                         *, smooth_obstructions=False, excluded_candidates=None):
    initialize()
    irr = list(tbl.Irr())
    one_dim = [(j+1,chi) for j,chi in enumerate(irr) if int(chi[0]) == 1]
    sizes = list(tbl.SizesConjugacyClasses())
    for parts in constituents([int(chi[0]) for chi in irr],n,reducible,max_representations):
        chi = sum((irr[i-1] for i in parts),0*irr[0])
        if int(libgap.GrassProjectiveKernelSize(tbl,chi)) != 1:
            continue
        sections = libgap.GrassSectionCharacter(tbl,chi,k,d)
        for j,psi in one_dim:
            multiplicity = int(libgap.ScalarProduct(tbl,sections,psi))
            if multiplicity <= 0:
                continue
            # A smooth connected section spans this representation. This is only
            # a conditional test; construct performs a pointwise restriction test.
            span_kernel = None
            if d == 1:
                span = sections-psi
                span_kernel = int(libgap.GrassProjectiveKernelSize(tbl,span))
            key = json.dumps([group_spec,k,n,d,parts,j],sort_keys=True)
            candidate = {
                'id':hashlib.sha256(key.encode()).hexdigest()[:12],
                'group':group_spec,'order':int(tbl.Size()),'k':k,'n':n,'d':d,
                'constituents':parts,'linear_character':j,
                'section_dimension':int(sections[0]),'multiplicity':multiplicity,
                'projective_kernel_order':1,
                'linear_span_projective_kernel_order':span_kernel,
                'status':'semi_invariant_section_exists; hypersurface_faithfulness_and_smoothness_unchecked',
                'character_values':[str(x) for x in chi],
                'linear_character_values':[str(x) for x in psi],
            }
            # The raw character catalogue remains available unchanged. Never
            # apply smooth-only geometry when singular sections are requested.
            if smooth_obstructions:
                obstruction = smooth_character_obstruction(tbl,chi,sections,psi,k,n,d)
                if obstruction is not None:
                    if excluded_candidates is not None:
                        excluded_candidates.append({'candidate':candidate,
                                                    'obstruction':obstruction})
                    continue
            yield candidate


def explicit_table(spec):
    initialize()
    if spec.startswith('small:'):
        _,order,index = spec.split(':')
        G=libgap.SmallGroup(int(order),int(index))
        return G.CharacterTable(),spec
    name = spec.removeprefix('table:')
    match = re.fullmatch(r'Alt\((\d+)\)', name)
    if match:
        name = 'A'+match[1]
    tbl=libgap.CharacterTable(name)
    if str(tbl)=='fail':
        raise ValueError(f'Unknown character table {name}')
    return tbl,'table:'+str(tbl.Identifier())


# Necessary conditions for SMOOTH complex hyperplane sections only.
# Benedetti-Manivel, Thm 1.1: every automorphism of a smooth Gr(3,10)
# hyperplane section extends to PGL(10): https://arxiv.org/abs/2110.00255
# Its projective stabilizer acts faithfully on the associated smooth
# Debarre-Voisin fourfold (Song, Lemma 3.1 and Prop. 6.1):
# https://arxiv.org/html/2106.13287v2
# A nonabelian simple group has no nontrivial character on H^{2,0}, so this
# action is symplectic. Hohn-Mason, Thm 4.9, gives these element orders:
# https://arxiv.org/html/1409.6055
# Consequently the only prime divisors of the group order are 2,3,5,7,11
# (Cauchy's theorem). Checking ALL class orders is stronger than that test.
# These conditions are NOT sufficient for smoothness, and are NOT a bound
# on group order. They do not justify truncating the catalogue at 100000.
DV_SYMPLECTIC_ORDERS = frozenset((1,2,3,4,5,6,7,8,9,11,12,14,15))


def smooth_group_obstruction(tbl, k, n, d):
    """Return a certified obstruction, or None (which is not a certificate).

    Gr(7,10) is included by Grassmannian duality. Cyclic simple groups are
    deliberately exempt: their action on H^{2,0} can be nontrivial. In
    particular neither 'order <= 22' nor the symplectic list applies to
    arbitrary cyclic actions. We do not infer a full-order cutoff from a
    totient bound on the order of the H^{2,0} character.
    """
    if n != 10 or k not in (3,7) or d != 1:
        return None
    if not bool(tbl.IsSimple()) or bool(tbl.IsAbelian()):
        return None
    forbidden = sorted({int(o) for o in tbl.OrdersClassRepresentatives()}
                       - DV_SYMPLECTIC_ORDERS)
    if not forbidden:
        return None
    return {'rule':'dv_nonabelian_simple_symplectic_orders',
            'forbidden_element_orders':forbidden,
            'allowed_element_orders':sorted(DV_SYMPLECTIC_ORDERS),
            'reason':'No smooth complex hyperplane section with a faithful action: '
                     'a nonabelian simple group must act symplectically on the '
                     'associated Debarre-Voisin fourfold.',
            'reference':'https://arxiv.org/html/1409.6055 (Theorem 4.9)'}


def search(args):
    validate(args.k,args.n,args.d)
    initialize()
    report={'schema':1,'parameters':vars(args).copy(),'scope':
            'Finite bounded GAP catalogue; ordinary linear representations only. Projective representations that do not lift are not enumerated.',
            'complete':False,'entries_examined':0,'tables_examined':0,'candidates':[],'errors':[],'excluded_groups':[],'excluded_candidates':[]}
    report['parameters'].pop('func',None)
    use_obstructions = getattr(args, 'smooth_obstructions', True)
    report['parameters']['smooth_obstructions'] = use_obstructions
    if use_obstructions:
        report['scope'] += (' Necessary smoothness obstructions enabled for '
                            'Gr(3,10)/Gr(7,10) hyperplanes: group orders, deformation '
                            'character positivity, and symplectic traces for perfect '
                            'groups. Excluded groups and candidates are recorded. '
                            'Surviving candidates are not smoothness certificates.')
    seen=set()
    if args.group:
        entries = args.group
    elif args.catalogue=='small':
        libgap.LoadPackage('smallgrp')
        def small_entries():
            for order in range(args.min_order,args.max_order+1):
                if not bool(libgap.SmallGroupsAvailable(order)):
                    report['errors'].append({'order':order,'error':'SmallGroups unavailable'})
                    continue
                for i in range(1,int(libgap.NumberSmallGroups(order))+1):
                    yield f'small:{order}:{i}'
        entries=small_entries()
    else:
        pred=libgap.eval(f'x -> x >= {int(args.min_order)} and x <= {int(args.max_order)}')
        if args.groups=='all':
            names=libgap.AllCharacterTableNames(libgap.Size,pred)
        else:
            names=libgap.AllCharacterTableNames(libgap.IsSimple,args.groups=='simple',libgap.Size,pred)
        entries=['table:'+str(x) for x in names]
    interrupted=False
    try:
        for spec in entries:
            if report['entries_examined']>=args.max_groups:
                interrupted=True
                break
            report['entries_examined']+=1
            try:
                tbl,spec=explicit_table(spec)
                if not args.group and args.catalogue=='tables' and bool(libgap.IsDuplicateTable(tbl)):
                    continue
                if spec in seen:
                    continue
                seen.add(spec)
                simple=bool(tbl.IsSimple())
                if args.groups=='simple' and not simple or args.groups=='nonsimple' and simple:
                    continue
                if not args.min_order <= int(tbl.Size()) <= args.max_order:
                    continue
                report['tables_examined']+=1
                if use_obstructions:
                    obstruction = smooth_group_obstruction(tbl,args.k,args.n,args.d)
                    if obstruction is not None:
                        report['excluded_groups'].append(dict(
                            group=spec, order=int(tbl.Size()), **obstruction))
                        if getattr(args,'verbose',False):
                            print(f'Excluded {spec}: {obstruction["reason"]}',flush=True)
                        continue
                if getattr(args,'verbose',False):
                    print(f'Searching {spec} (order {tbl.Size()})',flush=True)
                for candidate in candidates_for_table(
                        tbl,spec,args.k,args.n,args.d,args.reducible,args.max_representations,
                        smooth_obstructions=use_obstructions,
                        excluded_candidates=report['excluded_candidates']):
                    report['candidates'].append(candidate)
                    print(f"  {candidate['id']} ({spec}): constituents {candidate['constituents']}, character {candidate['linear_character']}, dimension {candidate['multiplicity']}",flush=True)
            except (RuntimeError,ValueError) as exc:
                report['errors'].append({'group':spec,'error':str(exc)})
                print(f'  Not completed: {exc}',flush=True)
            if args.output:
                write_json(args.output,report)
    except KeyboardInterrupt:
        interrupted=True
    report['complete']=not interrupted and not report['errors']
    report['stop_reason']='group limit or interruption' if interrupted else 'catalogue exhausted'
    if args.output:
        write_json(args.output,report)
    if args.output:
        print(f"Saved {len(report['candidates'])} candidates to {args.output}; complete within stated scope: {report['complete']}")
    return report
