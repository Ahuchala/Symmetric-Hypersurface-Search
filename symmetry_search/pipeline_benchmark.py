"""Compare complete smoothness invocations, including startup and checkpoints."""
import contextlib
from datetime import datetime, timezone
import hashlib
import io
import json
from pathlib import Path
import platform
import statistics
import subprocess
import tempfile
import time
from types import SimpleNamespace
from sage.env import SAGE_VERSION
from .common import ROOT, write_json
from .smoothness import smooth


def run_pipeline_benchmarks(*,model=None,prime=113,repeats=3,max_charts=12,seconds=60,output=None,verbose=False):
    """Fresh checkpoints per run. Set max_charts=None for the full A7 check.

    Compares dense/fresh reference, reusable worker, slimgb, contractions,
    and symmetry reduction. Contractions are included only for k=3,d=1.
    Partial runs must return incomplete, never smooth; full runs require smooth.
    Intended for the supplied smooth A7 fixture or another known smooth model.
    """
    path=Path(model) if model is not None else ROOT/'output/a7-equation.json'
    path=path.resolve();data=json.loads(path.read_text())
    if repeats<1 or int(repeats)!=repeats or seconds<=0 or (max_charts is not None and max_charts<1):raise ValueError('Invalid benchmark bounds')
    repeats=int(repeats);prime=int(prime);seconds=float(seconds)
    max_charts=None if max_charts is None else int(max_charts)
    variants=[('reference',dict(backend='fresh',dense=True)),
              ('persistent',dict(backend='persistent',dense=False)),
              ('parallel',dict(backend='persistent',dense=False,workers=2)),
              ('slimgb',dict(backend='persistent',dense=False,algorithm='slimgb'))]
    if (data['k'],data['d'])==(3,1):variants.append(('contraction',dict(backend='persistent',dense=False,algorithm='slimgb',method='contraction')))
    # Symmetry can cover more charts for the same representative budget.
    if max_charts is None:variants.append(('symmetry',dict(backend='persistent',dense=False,algorithm='slimgb',symmetry=True)))
    report={'schema':1,'created_utc':datetime.now(timezone.utc).isoformat(),
            'model_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
            'prime':prime,'repeats':repeats,'max_charts':max_charts,'timeout_seconds':seconds,
            'environment':{'sage':SAGE_VERSION,'python':platform.python_version(),'platform':platform.platform()},
            'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((ROOT/'symmetry_search').glob('*.py'))},
            'timing_policy':'Sequential variants; fresh checkpoint and worker pool per sample; wall time includes startup, proof verification, and checkpoint writes.',
            'variants':[],'passed':False}
    try:report['git_commit']=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    except (OSError,subprocess.CalledProcessError):report['git_commit']=None
    with tempfile.TemporaryDirectory(prefix='grass-pipeline-bench-') as tmp:
        for name,options in variants:
            row={'name':name,'options':options,'samples':[]}
            for rep in range(int(repeats)):
                args=dict(model=str(path),prime=prime,seconds=seconds,max_charts=1000000 if max_charts is None else max_charts,
                          output=str(Path(tmp)/f'{name}-{rep}.json'),restart=True,verbose=False,workers=1)
                args.update(options)
                start=time.perf_counter()
                with contextlib.redirect_stdout(io.StringIO()):evidence=smooth(SimpleNamespace(**args))
                elapsed=time.perf_counter()-start
                full=evidence['status']=='smooth_characteristic_zero'
                expected_count=min(max_charts,evidence['total_charts']) if max_charts is not None else evidence['total_charts']
                passed=(full if max_charts is None else
                        len(evidence['charts'])==expected_count and evidence['status']==('smooth_characteristic_zero' if expected_count==evidence['total_charts'] else 'incomplete'))
                passed=passed and all(c['status']=='smooth_chart' for c in evidence['charts'].values())
                row['samples'].append({'wall_seconds':elapsed,'passed':passed,'evidence':evidence})
                if verbose:print(name,rep+1,round(elapsed,3),'PASS' if passed else 'FAIL',flush=True)
            row['passed']=all(s['passed'] for s in row['samples'])
            row['median_wall_seconds']=statistics.median(s['wall_seconds'] for s in row['samples'])
            report['variants'].append(row)
            if output:write_json(output,report)
    report['passed']=all(r['passed'] for r in report['variants'])
    if report['passed']:
        reference=report['variants'][0]['median_wall_seconds']
        for row in report['variants']:row['speedup']=reference/row['median_wall_seconds']
    if output:write_json(output,report)
    return report
