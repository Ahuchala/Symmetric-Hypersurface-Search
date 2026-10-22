"""Exact chart checks, reusable killable workers, and certified chart orbits."""
import hashlib
import itertools
import json
import subprocess
import sys
import tempfile
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from .common import ROOT, write_json
from .models import read_model
from .chart_engine import ChartChecker, coefficient_reduction, check_chart, coordinate_singularity
from .chart_symmetry import chart_orbits
from .workers import PersistentWorker


def worker(args):
    model=read_model(args.model)
    try:
        result=check_chart(model,tuple(map(int,args.pivot.split(','))),args.prime,
            algorithm=getattr(args,'algorithm','auto'),method=getattr(args,'method','jacobian'),
            sparse=not getattr(args,'dense',False))
    except Exception as exc:result={'status':'error','error':str(exc)}
    write_json(args.output,result)


def fresh_check(path,pivot,prime,seconds,algorithm,method,sparse):
    start=time.perf_counter()
    with tempfile.TemporaryDirectory(prefix='grass-smooth-') as tmp:
        out=Path(tmp)/'result.json'
        command=[sys.executable,str(ROOT/'grassmannian_symmetry.py'),'_chart',str(path),
                 '--pivot',','.join(map(str,pivot)),'--prime',str(prime),'--output',str(out),
                 '--algorithm',algorithm,'--method',method]
        if not sparse:command.append('--dense')
        try:
            done=subprocess.run(command,capture_output=True,text=True,timeout=seconds)
            result=json.loads(out.read_text()) if done.returncode==0 and out.exists() else {'status':'error','error':done.stderr[-3000:]}
        except subprocess.TimeoutExpired:result={'status':'timeout','seconds':seconds}
    result['worker_wall_seconds']=time.perf_counter()-start
    return result


def smooth(args):
    started=time.perf_counter();path=Path(args.model).resolve();model=read_model(path)
    backend=getattr(args,'backend','fresh');workers=getattr(args,'workers',1)
    algorithm=getattr(args,'algorithm','auto');method=getattr(args,'method','jacobian')
    sparse=not getattr(args,'dense',False)
    if backend not in ('fresh','persistent'):raise ValueError('backend must be fresh or persistent')
    if int(workers)!=workers or workers<1:raise ValueError('workers must be a positive integer')
    workers=int(workers)
    if args.seconds<=0 or args.max_charts<0:raise ValueError('Invalid chart limit or timeout')
    checker=ChartChecker(model,args.prime,algorithm=algorithm,method=method,sparse=sparse) # validate before launch
    orbits,proofs=chart_orbits(model,args.prime,getattr(args,'symmetry',False),getattr(args,'symmetry_permutations',None))
    digest=hashlib.sha256(path.read_bytes()).hexdigest()
    report_path=Path(args.output or str(args.model)+f'.smooth-p{args.prime}.json')
    if report_path.exists() and not args.restart:
        report=json.loads(report_path.read_text())
        if report['model_sha256']!=digest or report['prime']!=args.prime:
            raise ValueError('Checkpoint belongs to another equation or prime; choose a new --output or --restart')
    else:
        report={'schema':1,'model':str(path),'model_sha256':digest,'prime':args.prime,'charts':{},'status':'incomplete'}
    pivots=list(itertools.combinations(range(model['n']),model['k']))
    report['total_charts']=len(pivots)
    report['chart_symmetry']={'generators':proofs,'orbits':[[list(I) for I in orbit] for orbit in orbits]}
    verbose=getattr(args,'verbose',False)
    passed=sum(c['status']=='smooth_chart' for c in report['charts'].values())
    field='characteristic zero' if args.prime==0 else f'F_{args.prime}'
    print(f'Smoothness check over {field}: {len(pivots)} charts, {len(orbits)} orbits; {passed} already passed.',flush=True)
    def key(I):return ','.join(map(str,I))
    def cover(orbit,representative):
        for I in orbit:
            if I!=representative and report['charts'].get(key(I),{}).get('status') not in ('smooth_chart','singular_chart','invalid_reduction'):
                report['charts'][key(I)]={'status':'smooth_chart','pivot':list(I),
                    'inferred_from':key(representative),'verification':'verified monomial chart symmetry',
                    'symmetry_generators':proofs}
    # Cheap global search for an exact singular coordinate-plane witness.
    # A zero chart budget remains a request to inspect checkpoints only.
    if args.max_charts>0 and not any(c['status'] in ('singular_chart','invalid_reduction') for c in report['charts'].values()) and passed<len(pivots):
        witness=coordinate_singularity(checker)
        if witness is not None:
            report['charts'][key(witness['pivot'])]=witness
            report['coordinate_singularity_precheck']=witness
            print('Exact singular coordinate plane:',tuple(witness['pivot']),flush=True)
    pending=[]
    for orbit in orbits:
        known=next((I for I in orbit if report['charts'].get(key(I),{}).get('status')=='smooth_chart'),None)
        if known is not None:cover(orbit,known)
        else:pending.append(orbit)
    # A previously certified singular/zero chart already settles the check.
    if any(c['status'] in ('singular_chart','invalid_reduction') for c in report['charts'].values()):pending=[]
    pending=pending[:args.max_charts]
    pool=[PersistentWorker(path,args.prime,algorithm=algorithm,method=method,sparse=sparse)
          for _ in range(min(workers,len(pending)))] if backend=='persistent' else []
    attempted=0
    try:
        with ThreadPoolExecutor(max_workers=workers) as executor:
            for offset in range(0,len(pending),workers):
                batch=pending[offset:offset+workers]
                if backend=='persistent':
                    futures=[executor.submit(pool[i].check,orbit[0],args.seconds) for i,orbit in enumerate(batch)]
                else:
                    futures=[executor.submit(fresh_check,path,orbit[0],args.prime,args.seconds,algorithm,method,sparse) for orbit in batch]
                stop=False
                for orbit,future in zip(batch,futures):
                    result=future.result();attempted+=1
                    report['charts'][key(orbit[0])]=result
                    if result['status']=='smooth_chart':cover(orbit,orbit[0])
                    if result['status'] in ('singular_chart','invalid_reduction'):stop=True
                    write_json(report_path,report)
                    if verbose:print(orbit[0],result['status'],flush=True)
                    elif attempted%10==0:print(f'  {attempted} chart representatives checked.',flush=True)
                if stop:break
    finally:
        for process in pool:process.close()
    report['last_run']={'backend':backend,'workers':workers,'algorithm':algorithm,'method':method,
                        'sparse_minors':sparse,'charts_computed':attempted,
                        'worker_starts':sum(p.starts for p in pool) if pool else attempted,
                        'wall_seconds':time.perf_counter()-started}
    statuses=[x['status'] for x in report['charts'].values()]
    if 'invalid_reduction' in statuses:
        report['status']='invalid_reduction'
        report['characteristic_zero_conclusion']='inconclusive'
    elif 'singular_chart' in statuses:
        report['status']='singular_characteristic_zero' if args.prime==0 else 'singular_special_fibre'
        report['characteristic_zero_conclusion']='singular' if args.prime==0 else 'inconclusive; a singular reduction does not prove characteristic-zero singularity'
    elif len(statuses)==len(pivots) and all(s=='smooth_chart' for s in statuses):
        report['status']='smooth_characteristic_zero'
        report['characteristic_zero_conclusion']='smooth'
        report['proof']='All affine Grassmannian charts have unit Jacobian ideal'+('; the nonzero integral section defines a flat projective Cartier divisor over the local coefficient ring, so smooth special fibre implies smooth generic fibre' if args.prime else '')
    else:
        report['status']='incomplete'
        report['characteristic_zero_conclusion']='inconclusive'
    write_json(report_path,report)
    if verbose:
        print(report['status']+': '+report['characteristic_zero_conclusion'])
    elif report['status']=='smooth_characteristic_zero':
        print(f'Smooth in characteristic zero. All {len(pivots)}/{len(pivots)} charts passed.')
    elif report['status']=='singular_characteristic_zero':
        print('Singular in characteristic zero.')
    elif report['status']=='singular_special_fibre':
        print(f'Singular over F_{args.prime}; characteristic-zero smoothness is inconclusive.')
    elif report['status']=='invalid_reduction':
        print('Invalid reduction: the equation reduces to zero. Choose another prime.')
    else:
        counts=Counter(statuses)
        issues=', '.join(f"{counts[status]} {label}" for status,label in [('timeout','timed out'),('error','errors')] if counts[status])
        detail=f'; {issues}' if issues else ''
        print(f"Incomplete: {counts['smooth_chart']}/{len(pivots)} charts passed{detail}. Rerun to resume.")
    print('Evidence:',report_path)
    return report
