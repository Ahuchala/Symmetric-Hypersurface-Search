#!/usr/bin/env sage -python
"""Search -> construct sparse semi-invariant equations -> verify smoothness.
Run with `sage -python grassmannian_symmetry.py --help`.
"""
import argparse
import json
import sys
import os
import subprocess
from pathlib import Path


def parser():
    p=argparse.ArgumentParser(description=__doc__)
    sub=p.add_subparsers(dest='command',required=True)
    s=sub.add_parser('search',help='Find projectively faithful representations with nonzero semi-invariant Grassmannian sections')
    for x in ('k','n','d'):s.add_argument(x,type=int)
    s.add_argument('--groups',choices=['simple','nonsimple','all'],default='simple')
    s.add_argument('--reducible',action='store_true',help='Include all degree-n sums of irreducibles within the stated limit')
    s.add_argument('--catalogue',choices=['tables','small'],default='tables')
    s.add_argument('--group',action='append',help='Restrict to a table name (e.g. A7, L2(11)) or small:order:index; repeatable')
    s.add_argument('--min-order',type=int,default=2)
    s.add_argument('--max-order',type=int,default=100000)
    s.add_argument('--max-groups',type=int,default=1000)
    s.add_argument('--max-representations',type=int,default=10000)
    s.add_argument('--no-smooth-obstructions',dest='smooth_obstructions',action='store_false',
                   help='Retain groups obstructed for smooth DV hyperplanes (e.g. to study singular sections)')
    s.add_argument('--output',default='output/symmetry-search.json')
    c=sub.add_parser('construct',help='Select a search candidate and construct an exactly verified sparse equation')
    c.add_argument('search');c.add_argument('candidate')
    c.add_argument('--output',default='output/symmetry-equation.json')
    c.add_argument('--basis-tries',type=int,default=2,help='Additional cyclic eigenbases to try; 0 keeps original and known certified bases')
    c.add_argument('--max-conductor',type=int,default=420)
    c.add_argument('--max-sections',type=int,default=None,help='Explicit section bound; default 1000 for expanded construction, entry-budget bounds for quadrics')
    c.add_argument('--quadratic-method',choices=['auto','reynolds','expanded'],default='auto')
    c.add_argument('--quadratic-trials',type=int,default=0,help='Bound on Reynolds seeds; 0 uses a dimension-dependent bound')
    c.add_argument('--max-entries',type=int,default=2000000)
    c.add_argument('--trials',type=int,default=20,help='Small combinations of invariant basis vectors to test for sparsity')
    c.add_argument('--seed',type=int,default=0)
    c.add_argument('--seconds',type=float,default=300,help='Hard time limit for the complete construction; 0 disables')
    c.add_argument('--coefficients',help='Comma-separated coefficients in the original-basis invariant space')
    c.add_argument('--random-section',action='store_true',help='Use a reproducible small random combination instead of preferring sparse basis vectors')
    c.add_argument('--gap-group',help='Explicit trusted GAP group expression for a table without an available group model')
    c.add_argument('--faithfulness-charts',type=int,default=1,help='Charts to use for exact pointwise faithfulness; 0 leaves unchecked')
    f=sub.add_parser('faithful',help='Check faithfulness on the reduced hypersurface itself')
    f.add_argument('model');f.add_argument('--max-charts',type=int,default=1000000)
    f.add_argument('--output',help='Evidence file; defaults to MODEL.faithful.json')
    s=sub.add_parser('smooth',help='Check all Grassmannian charts, with resumable bounded workers')
    s.add_argument('model');s.add_argument('--prime',type=int,default=0,help='0 for characteristic zero; a good prime can certify smoothness, but not characteristic-zero singularity')
    s.add_argument('--seconds',type=float,default=120,help='Timeout per chart, in seconds')
    s.add_argument('--max-charts',type=int,default=1000000,help='Maximum new charts this invocation; partial coverage is never reported as smooth')
    s.add_argument('--output');s.add_argument('--restart',action='store_true')
    show=sub.add_parser('show',help='Print the saved equation and its status');show.add_argument('model');show.add_argument('--latex',action='store_true')
    w=sub.add_parser('_chart',help=argparse.SUPPRESS)
    w.add_argument('model');w.add_argument('--pivot',required=True);w.add_argument('--prime',type=int,required=True);w.add_argument('--output',required=True)
    for name in ('smooth','_chart'):
        sub.choices[name].add_argument('--algorithm',choices=['auto','std','slimgb'],default='auto')
        sub.choices[name].add_argument('--method',choices=['jacobian','contraction'],default='jacobian')
        sub.choices[name].add_argument('--dense',action='store_true',help='Compute all minors (reference path)')
    sub.choices['smooth'].add_argument('--backend',choices=['fresh','persistent'],default='persistent')
    sub.choices['smooth'].add_argument('--workers',type=int,default=1)
    sub.choices['smooth'].add_argument('--symmetry',action='store_true',help='Verify monomial symmetries and reduce chart orbits')
    for command in ('search','construct','faithful','smooth'):
        sub.choices[command].add_argument('--verbose',action='store_true',help='Show detailed progress and diagnostics')
    return p


def main():
    args=parser().parse_args()
    try:
        if args.command=='construct' and args.seconds>0 and os.environ.get('GRASS_SYMMETRY_CONSTRUCT_CHILD')!='1':
            env=os.environ.copy();env['GRASS_SYMMETRY_CONSTRUCT_CHILD']='1'
            try:
                return subprocess.run([sys.executable,str(Path(__file__).resolve()),*sys.argv[1:]],env=env,timeout=args.seconds).returncode
            except subprocess.TimeoutExpired:
                print('Construction time limit reached. No unfinished stage is certified; increase --seconds to retry.',file=sys.stderr)
                return 2
        if args.command=='search':
            from symmetry_search.search import search
            search(args)
        elif args.command=='construct':
            from symmetry_search.construct import construct
            construct(args)
        elif args.command in ('smooth','_chart'):
            from symmetry_search.smoothness import smooth,worker
            (worker if args.command=='_chart' else smooth)(args)
        elif args.command=='show':
            from symmetry_search.models import read_model,print_equation
            model=read_model(args.model)
            if args.latex:
                from symmetry_search.models import latex_equation
                print(latex_equation(model))
            else:
                print_equation(model)
        elif args.command=='faithful':
            from symmetry_search.models import read_model
            from symmetry_search.faithfulness import check_faithfulness
            from symmetry_search.common import write_json
            result=check_faithfulness(read_model(args.model),args.max_charts)
            evidence_path=args.output or args.model+'.faithful.json'
            write_json(evidence_path,result)
            if args.verbose:
                print(json.dumps(result,indent=2))
            else:
                print('Faithfulness:',result['status'])
                print('Evidence:',evidence_path)
    except (ValueError,RuntimeError,ArithmeticError) as exc:
        print('Not completed:',exc,file=sys.stderr)
        return 2
    return 0

if __name__=='__main__':
    sys.exit(main())
