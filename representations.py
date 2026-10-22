"""Construct selected ordinary representations over explicit cyclotomic fields."""
import json
import re
from sage.all import libgap, matrix, block_diagonal_matrix, sage_eval, QQ, CyclotomicField
from .common import ROOT, field_from
from .search import initialize, explicit_table


def existing_matrices(name):
    if name=='A7':
        K=CyclotomicField(7,names='z');z=K.gen();b=z+z**2+z**4;B=-1-b
        path=ROOT/'A7Intersection.M2';context={'b':b,'B':B}
    elif name in ('L2(11)','PSL211'):
        K=QQ;path=ROOT/'PSL211Intersection.M2';context={}
    else:
        return None
    text=path.read_text()
    mats=[]
    for letter in ('r','s'):
        expr=re.search(letter+r' = matrix\(K,\s*(\{.*?\})\);',text,re.S).group(1)
        mats.append(matrix(K,sage_eval(expr.replace('{','[').replace('}',']'),locals=context)))
    return K,mats


def group_for(spec):
    if spec.startswith('small:'):
        _,order,index=spec.split(':')
        return libgap.SmallGroup(int(order),int(index)),None
    name=spec.removeprefix('table:')
    builtin=existing_matrices(name)
    if builtin:
        K,mats=builtin
        return libgap.Group([libgap(M) for M in mats]),builtin
    if re.fullmatch(r'A\d+',name):
        return libgap.AlternatingGroup(int(name[1:])),None
    if re.fullmatch(r'S\d+',name):
        return libgap.SymmetricGroup(int(name[1:])),None
    match=re.fullmatch(r'L(\d+)\((\d+)\)',name)
    if match:
        return libgap.PSL(int(match[1]),int(match[2])),None
    match=re.fullmatch(r'U(\d+)\((\d+)\)',name)
    if match:
        return libgap.PSU(int(match[1]),int(match[2])),None
    match=re.fullmatch(r'M(11|12|22|23|24)',name)
    if match:
        return libgap.MathieuGroup(int(match[1])),None
    libgap.LoadPackage('atlasrep')
    G=libgap.AtlasGroup(name)
    if str(G)=='fail':
        raise ValueError(f'No explicit group model available for {spec}. Use the SmallGroups catalogue or provide --gap-group for this table.')
    return G,None


def realize(candidate,gap_group=None):
    initialize()
    source_table,_=explicit_table(candidate['group'])
    source_irreps=list(source_table.Irr())
    source_chi=sum((source_irreps[i-1] for i in candidate['constituents']),0*source_irreps[0])
    source_psi=source_irreps[candidate['linear_character']-1]
    if [str(x) for x in source_chi]!=candidate['character_values'] or [str(x) for x in source_psi]!=candidate['linear_character_values']:
        raise ValueError('Character-table ordering or values changed since the saved search; rerun the search before selecting this candidate')
    if gap_group:
        G=libgap.eval(gap_group);builtin=None
        if not bool(libgap.IsGroup(G)):
            raise ValueError('--gap-group must evaluate to a finite GAP group')
    else:
        G,builtin=group_for(candidate['group'])
    concrete=G.CharacterTable()
    change=libgap.TransformingPermutationsCharacterTables(source_table,concrete)
    if str(change)=='fail':
        raise ValueError('Constructed group does not have the selected character table')
    rowperm=change['rows']
    indices=[int(libgap.OnPoints(i,rowperm))-1 for i in candidate['constituents']]
    linear_index=int(libgap.OnPoints(candidate['linear_character'],rowperm))-1
    irreps=list(concrete.Irr())
    target=sum((irreps[i] for i in indices),0*irreps[0]);psi=irreps[linear_index]
    classes=list(G.ConjugacyClasses())
    reps=[cl.Representative() for cl in classes]
    gens=list(G.GeneratorsOfGroup())
    psi_gens=[libgap.GrassCharacterValue(G,psi,g) for g in gens]
    psi_classes=list(psi)
    rawgen=None;rawclass=None;seed_bases=[]
    if builtin and candidate['n']==10:
        K,natural=builtin
        for conjugate in (False,True):
            converter=lambda x: K(x.sage()).conjugate() if conjugate else K(x.sage())
            classm=[matrix(K,[[converter(x) for x in row] for row in g]) for g in reps]
            char=[libgap(M.trace()) for M in classm]
            if all(x==y for x,y in zip(char,target)):
                mats=[matrix(K,[[converter(x) for x in row] for row in g]) for g in gens]
                if candidate['group']=='table:A7':
                    certificate=json.loads((ROOT/'symmetry_search/data/a7-sparse-certificate.json').read_text())
                    C=matrix(K,[[sage_eval(x,locals={'zeta7':K.gen()}) for x in row] for row in certificate['change_of_basis']])
                    if conjugate:
                        C=C.apply_map(lambda x:x.conjugate())
                    seed_bases.append(('certified-A7-cyclic-basis',C))
                break
        else:
            mats=None
    else:
        mats=None
    if mats is None:
        rho=[libgap.IrreducibleRepresentationsDixon(G,irreps[i]) for i in sorted(set(indices))]
        lookup=dict(zip(sorted(set(indices)),rho))
        rawgen=[[list(map(list,libgap.Image(lookup[i],g))) for i in indices] for g in gens]
        rawclass=[[list(map(list,libgap.Image(lookup[i],g))) for i in indices] for g in reps]
        scalars=[x for images in rawgen+rawclass for block in images for row in block for x in row]+psi_gens+psi_classes
        K=field_from(int(libgap.Conductor(scalars)))
        def convert(images):
            return block_diagonal_matrix([matrix(K,[[K(x.sage()) for x in row] for row in block]) for block in images])
        mats=[convert(m) for m in rawgen];classm=[convert(m) for m in rawclass]
    pg=[K(x.sage()) for x in psi_gens];pc=[K(x.sage()) for x in psi_classes]
    assert all(M.nrows()==candidate['n'] for M in mats)
    assert all(libgap(M.trace())==v for M,v in zip(classm,target))
    return {'field':K,'generators':mats,'generator_characters':pg,
            'classes':[{'matrix':M,'character':p,'size':int(cl.Size()),'order':int(g.Order())} for M,p,cl,g in zip(classm,pc,classes,reps)],
            'seed_bases':seed_bases,'group_order':int(G.Size())}
