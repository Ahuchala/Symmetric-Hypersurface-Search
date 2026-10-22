"""Fresh-process correctness/timing benchmarks; callable from a Sage notebook."""
import argparse
from datetime import datetime, timezone
import hashlib
import itertools
import json
import platform
from pathlib import Path
import statistics
import subprocess
import sys
import tempfile
import time

from sage.env import SAGE_VERSION
from .common import ROOT, write_json
from .smoothness import check_chart


def model(k,n,d,terms):
    """Small explicit fixture, independent of character search/construction."""
    return {'schema':1,'k':k,'n':n,'d':d,'field_conductor':1,
            'plucker_indices':[list(I) for I in itertools.combinations(range(n),k)],
            'terms':[{'monomial':list(m),'coefficient':[str(c)]} for m,c in terms]}


def fingerprint(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def cases(profile='quick'):
    """Quick includes three A7 charts; full includes all 120. No cached results."""
    if profile not in ('quick','full'):raise ValueError('profile must be quick or full')
    fixtures=[
        ('gr24-symplectic',model(2,4,1,[((0,),1),((5,),1)]),0,
         list(itertools.combinations(range(4),2)),'smooth_chart'),
        # The singular point is span(e2,e3); it is missed by the first big cell.
        ('gr24-boundary-singular',model(2,4,1,[((0,),1)]),0,[(2,3)],'singular_chart'),
        ('conic-Q',model(1,3,2,[((0,0),1),((1,1),1),((2,2),1)]),0,[(0,)],'smooth_chart'),
        ('conic-bad-prime',model(1,3,2,[((0,0),1),((1,1),1),((2,2),1)]),2,[(0,)],'singular_chart'),
        ('nonreduced-point',model(1,2,2,[((0,0),1)]),0,[(1,)],'singular_chart'),
        ('empty-chart',model(1,2,1,[((0,),1)]),0,[(0,)],'smooth_chart'),
        ('zero-reduction',model(1,2,1,[((0,),3)]),3,[(0,)],'invalid_reduction'),
    ]
    a7=json.loads((ROOT/'output/a7-equation.json').read_text())
    pivots=list(itertools.combinations(range(10),3)) if profile=='full' else [(0,1,2),(1,4,7),(7,8,9)]
    fixtures.append(('a7-p113',a7,113,pivots,'smooth_chart'))
    result=[]
    for name,equation,prime,pivots,expected in fixtures:
        for pivot in pivots:
            result.append({'id':name+':'+','.join(map(str,pivot)), 'model':equation,
                           'prime':prime,'pivot':list(pivot),'expected':expected})
    return result


def _worker(source,target):
    case=json.loads(Path(source).read_text())
    start=time.perf_counter()
    result=check_chart(case['model'],tuple(case['pivot']),case['prime'])
    elapsed=time.perf_counter()-start
    write_json(target,{'result':result,'calculation_seconds':elapsed})


def compare_reports(current,baseline):
    """Reject changed fixtures/settings or failed runs before comparing medians."""
    if not current['passed'] or not baseline['passed']:
        raise ValueError('Both benchmark reports must pass correctness checks')
    for key in ('schema','suite_fingerprint','settings'):
        if current[key]!=baseline[key]:raise ValueError(f'Incompatible benchmark {key}')
    for key in ('sage','python','system','machine'):
        if current['environment'][key]!=baseline['environment'][key]:
            raise ValueError(f'Incompatible environment: {key}')
    old={r['id']:r for r in baseline['cases']}
    if set(old)!={r['id'] for r in current['cases']}:raise ValueError('Different case sets')
    return [{'id':r['id'],
             'calculation_speedup':old[r['id']]['median_calculation_seconds']/r['median_calculation_seconds'],
             'wall_speedup':old[r['id']]['median_wall_seconds']/r['median_wall_seconds']}
            for r in current['cases']]


def run_benchmarks(*,profile='quick',repeats=3,seconds=60,output=None,baseline=None,verbose=False):
    """Time fresh bounded workers; return evidence, medians, and correctness status.

    Quick is NOT an all-chart smoothness certificate for A7. Full recomputes
    every A7 chart. Wall time includes Sage startup; calculation time excludes it.
    Failed/time-limited cases remain failures, not speed improvements.
    """
    if int(repeats)!=repeats or repeats<1 or seconds<=0:raise ValueError('Invalid repeats or seconds')
    if isinstance(baseline,(str,Path)) and output and Path(baseline).resolve()==Path(output).resolve():
        raise ValueError('Use a different output path to preserve the baseline')
    selected=cases(profile)
    try:revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    except (OSError,subprocess.CalledProcessError):revision=None
    sources={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
             for p in sorted((ROOT/'symmetry_search').glob('*.py'))}
    report={'schema':1,'created_utc':datetime.now(timezone.utc).isoformat(),
            'settings':{'profile':profile,'repeats':int(repeats),'timeout_seconds':float(seconds)},
            'suite_fingerprint':fingerprint(selected),'git_commit':revision,'source_sha256':sources,
            'environment':{'sage':SAGE_VERSION,'python':platform.python_version(),
                           'system':platform.platform(),'machine':platform.machine(),'processor':platform.processor()},
            'timing_policy':'Fresh worker each sample; no checkpoint reuse; no warmup; sequential cases.',
            'passed':False,'cases':[]}
    with tempfile.TemporaryDirectory(prefix='symmetry-bench-') as tmp:
        source=Path(tmp)/'case.json';target=Path(tmp)/'result.json'
        for case in selected:
            row={k:case[k] for k in ('id','prime','pivot','expected')}
            row['model_sha256']=fingerprint(case['model']);row['samples']=[]
            write_json(source,case)
            for _ in range(int(repeats)):
                target.unlink(missing_ok=True)
                start=time.perf_counter()
                try:
                    done=subprocess.run([sys.executable,'-m','symmetry_search.benchmark','--worker',str(source),str(target)],
                                        cwd=ROOT,capture_output=True,text=True,timeout=seconds)
                    wall=time.perf_counter()-start
                    if done.returncode or not target.exists():
                        sample={'status':'error','error':done.stderr[-3000:],'wall_seconds':wall}
                    else:
                        payload=json.loads(target.read_text())
                        sample={'status':payload['result']['status'],'evidence':payload['result'],
                                'calculation_seconds':payload['calculation_seconds'],'wall_seconds':wall}
                except subprocess.TimeoutExpired:
                    sample={'status':'timeout','wall_seconds':time.perf_counter()-start}
                row['samples'].append(sample)
            row['passed']=all(s['status']==case['expected'] for s in row['samples'])
            if row['passed']:
                row['median_calculation_seconds']=statistics.median(s['calculation_seconds'] for s in row['samples'])
                row['median_wall_seconds']=statistics.median(s['wall_seconds'] for s in row['samples'])
            report['cases'].append(row)
            if output:write_json(output,report)
            if verbose:print(case['id'], 'PASS' if row['passed'] else 'FAIL',flush=True)
    report['passed']=all(r['passed'] for r in report['cases'])
    if baseline:
        previous=json.loads(Path(baseline).read_text()) if isinstance(baseline,(str,Path)) else baseline
        report['comparison']=compare_reports(report,previous)
    if output:write_json(output,report)
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--profile',choices=['quick','full'],default='quick')
    parser.add_argument('--repeats',type=int,default=3)
    parser.add_argument('--seconds',type=float,default=60)
    parser.add_argument('--output',default='output/symmetry-benchmark.json')
    parser.add_argument('--baseline')
    parser.add_argument('--worker',nargs=2,help=argparse.SUPPRESS)
    args=parser.parse_args()
    if args.worker:_worker(*args.worker);return
    report=run_benchmarks(profile=args.profile,repeats=args.repeats,seconds=args.seconds,
                          output=args.output,baseline=args.baseline,verbose=True)
    print(('PASS' if report['passed'] else 'FAIL')+f": {len(report['cases'])} cases; evidence: {args.output}")
    if not report['passed']:raise SystemExit(1)

if __name__=='__main__':main()
