"""Interactive Sage interface; the command-line interface uses the same engines."""
from collections.abc import Sequence
from contextlib import contextmanager
from copy import deepcopy
from functools import cached_property
from html import escape
from pathlib import Path
from types import SimpleNamespace
import contextlib
import json
import os
import tempfile

from sage.all import PolynomialRing, prod
from .common import decode, decode_matrix, field_from, write_json
from .search import search as _run_search


@contextmanager
def _messages(verbose):
    if verbose:
        yield
    else:
        # Do not accumulate a potentially long search's output in notebook memory.
        with open(os.devnull,'w') as sink, contextlib.redirect_stdout(sink):
            yield


def search(k,n,d,*,group=None,simple=True,reducible=False,catalogue='tables',
           min_order=2,max_order=100000,max_groups=1000,max_representations=10000,
           verbose=False,smooth_obstructions=True):
    """Search a bounded catalogue, returning indexable SearchResults.

    simple=True: simple groups; False: non-simple groups; None: both.
    group may be a name, e.g. 'A7', or a list of names. Irreducibles are the
    default; reducible=True also includes sums. Hits are candidates, not
    smoothness or hypersurface-faithfulness certificates.
    smooth_obstructions=True applies proven necessary conditions for smooth
    Gr(3,10)/Gr(7,10) hyperplanes. Rejections are in report['excluded_groups']
    and report['excluded_candidates']. Set False to disable ALL smooth-only
    pruning and retain candidates for possibly singular sections. Cyclic
    groups are exempt from symplectic restrictions, not character positivity.
    """
    if simple is not None and type(simple) is not bool:
        raise ValueError('simple must be True, False, or None')
    if type(smooth_obstructions) is not bool:
        raise ValueError('smooth_obstructions must be True or False')
    if catalogue not in ('tables','small'):
        raise ValueError("catalogue must be 'tables' or 'small'")
    groups=None if group is None else ([group] if isinstance(group,str) else list(group))
    if groups is not None and not groups:
        raise ValueError('group cannot be an empty list')
    args=SimpleNamespace(k=int(k),n=int(n),d=int(d),group=groups,
        groups='all' if simple is None else ('simple' if simple else 'nonsimple'),
        reducible=bool(reducible),catalogue=catalogue,min_order=int(min_order),
        max_order=int(max_order),max_groups=int(max_groups),
        max_representations=int(max_representations),output=None,verbose=verbose,
        smooth_obstructions=smooth_obstructions)
    if args.min_order<1 or args.max_order<args.min_order or args.max_groups<1 or args.max_representations<1:
        raise ValueError('Search bounds must be positive and min_order <= max_order')
    with _messages(verbose):
        report=_run_search(args)
    return SearchResults(report)


class SearchResults(Sequence):
    """A list-like search result with explicit coverage and error information."""
    def __init__(self,report):
        self._report=deepcopy(report)
        self._candidates=tuple(Candidate(c) for c in report['candidates'])
        self.last_smooth_search=None

    def __len__(self):return len(self._candidates)

    def __getitem__(self,index):
        if isinstance(index,str):
            matches=[c for c in self if c.id==index]
            if not matches:raise KeyError(index)
            return matches[0]
        return self._candidates[index]

    @property
    def complete(self):return self._report['complete']

    @property
    def errors(self):return deepcopy(self._report.get('errors',[]))

    @property
    def report(self):return deepcopy(self._report)

    def __repr__(self):
        scope='complete within bounds' if self.complete else 'INCOMPLETE; inspect .errors and .report'
        lines=[f'{len(self)} candidates ({scope})']
        excluded=self._report.get('excluded_groups',[])
        if excluded:
            lines.append(f'{len(excluded)} groups excluded by proven smoothness '
                         'obstructions; see .report["excluded_groups"]')
        excluded_candidates=self._report.get('excluded_candidates',[])
        if excluded_candidates:
            lines.append(f'{len(excluded_candidates)} candidate families excluded by '
                         'proven smoothness obstructions; '
                         'see .report["excluded_candidates"]')
        for i,c in enumerate(self._candidates[:20]):
            lines.append(f'[{i}] {c}')
        if len(self)>20:lines.append(f'... {len(self)-20} more; access by index or iterate')
        return '\n'.join(lines)

    def _repr_html_(self):return '<pre>'+escape(repr(self))+'</pre>'

    def find_smooth(self, *, sections_per_candidate=5, max_candidates=None,
                    seed=0, prime=113, seconds=10, max_charts=1000000,
                    algorithm="auto", symmetry=False, workers=1,
                    faithfulness_charts=1000000, construction_options=None,
                    verbose=False):
        """Return the first certified smooth, faithful Hypersurface, or None.

        Try candidates in order, first the sparse equation, then reproducible
        combinations of its cached invariant basis. One-dimensional spaces
        need only one equation. Limits/timeouts are inconclusive, not rejection
        of the family. Inspect last_smooth_search for every attempted equation.
        seconds is per chart; construction has no wall-clock timeout.
        prime=0 checks characteristic zero directly. A smooth good reduction
        also certifies characteristic-zero smoothness.
        """
        from .smooth_search import find_smooth
        return find_smooth(self, sections_per_candidate=sections_per_candidate,
            max_candidates=max_candidates, seed=seed, prime=prime, seconds=seconds,
            max_charts=max_charts, algorithm=algorithm, symmetry=symmetry,
            workers=workers, faithfulness_charts=faithfulness_charts,
            construction_options=construction_options, verbose=verbose)

    def save(self,path):
        write_json(path,self._report)
        return Path(path)


class Candidate:
    """One group, representation, and scalar character from a search."""
    def __init__(self,data):self._data=deepcopy(data)

    @property
    def id(self):return self._data['id']

    @property
    def group(self):return self._data['group']

    @property
    def multiplicity(self):return self._data['multiplicity']

    @property
    def data(self):return deepcopy(self._data)

    def __repr__(self):
        c=self._data
        return f"{self.group}: representation {c['constituents']}, character {c['linear_character']}, invariant space dimension {self.multiplicity}"

    def construct(self,*,basis_tries=2,max_conductor=420,max_sections=None,
                  max_entries=2000000,trials=20,seed=0,coefficients=None,
                  random_section=False,gap_group=None,faithfulness_charts=1,
                  verbose=False,quadratic_method="auto",quadratic_trials=0):
        """Construct in memory, returning a Hypersurface with a Sage polynomial.

        Construction runs in the calling Sage process, with dimension/entry
        bounds (no hard wall-clock limit). Interrupt the kernel to stop it.
        coefficients is a list in the original-basis invariant space.
        Quadrics use Reynolds forms and exact Plucker straightening by default.
        quadratic_method="expanded" restores the older invariance system.
        max_sections=None uses 1000 for the expanded method; quadrics instead
        use max_entries for orbit, invariant-basis, and straightening bounds.
        Cyclic basis_tries/max_conductor apply to the expanded method only.
        """
        from .construct import construct_candidate
        if coefficients is not None and random_section:
            raise ValueError('Choose coefficients or random_section, not both')
        if coefficients is not None:
            coefficients=','.join(map(str,coefficients))
        args=SimpleNamespace(basis_tries=int(basis_tries),max_conductor=int(max_conductor),
            max_sections=None if max_sections is None else int(max_sections),max_entries=int(max_entries),trials=int(trials),
            seed=int(seed),coefficients=coefficients,random_section=bool(random_section),
            gap_group=gap_group,faithfulness_charts=int(faithfulness_charts),
            output=None,verbose=verbose,quadratic_method=quadratic_method,quadratic_trials=int(quadratic_trials))
        with _messages(verbose):
            model=construct_candidate(self._data,args)
        return Hypersurface(model)


class CheckResult:
    """A check's conclusion and evidence; unknown is never coerced to True."""
    def __init__(self,report,kind):
        self._report=deepcopy(report);self.kind=kind

    @property
    def status(self):return self._report['status']

    @property
    def report(self):return deepcopy(self._report)

    @property
    def is_smooth(self):
        if self.kind!='smoothness':raise AttributeError('This is not a smoothness check')
        return True if self.status=='smooth_characteristic_zero' else False if self.status=='singular_characteristic_zero' else None

    @property
    def is_faithful(self):
        if self.kind!='faithfulness':raise AttributeError('This is not a faithfulness check')
        return True if self.status=='faithful' else False if self.status=='not_faithful' else None

    def __bool__(self):
        raise TypeError('Use .is_smooth or .is_faithful explicitly: None means inconclusive')

    def __repr__(self):
        if self.kind=='smoothness':
            charts=self._report.get('charts',{})
            passed=sum(c['status']=='smooth_chart' for c in charts.values())
            total=self._report.get('total_charts','?')
            return f'{self.status} ({passed}/{total} charts passed)'
        return f'Faithfulness: {self.status}'

    def _repr_html_(self):return '<pre>'+escape(repr(self))+'</pre>'


class Hypersurface:
    """A verified equation with Sage polynomial, matrices, and optional checks."""
    def __init__(self,model,path=None):
        self._model=deepcopy(model)
        self._path=Path(path).resolve() if path is not None else None
        self._scratch=None
        self.last_smoothness=None
        self.last_faithfulness=None

    @property
    def k(self):return self._model['k']
    @property
    def n(self):return self._model['n']
    @property
    def d(self):return self._model['d']
    @property
    def data(self):return deepcopy(self._model)

    @cached_property
    def field(self):return field_from(self._model['field_conductor'])

    @cached_property
    def ring(self):
        names=['p_'+'_'.join(map(str,I)) for I in self._model['plucker_indices']]
        return PolynomialRing(self.field,names,order='degrevlex')

    def _polynomial(self,terms):
        return sum((decode(t['coefficient'],self.field)*prod(self.ring.gen(j) for j in t['monomial']) for t in terms),self.ring.zero())

    @cached_property
    def polynomial(self):return self._polynomial(self._model['terms'])

    @property
    def f(self):return self.polynomial

    @cached_property
    def invariant_basis(self):return tuple(self._polynomial(t) for t in self._model['invariant_basis'])

    @property
    def generators(self):return tuple(decode_matrix(M,self.field) for M in self._model['generators'])

    @property
    def change_of_basis(self):return decode_matrix(self._model['change_of_basis'],self.field)

    def p(self,*indices):
        """Return an alternating Plucker coordinate, e.g. X.p(0,1,7)."""
        if len(indices)==1 and isinstance(indices[0],(tuple,list)):
            indices=tuple(indices[0])
        if len(indices)!=self.k or any(int(i)!=i or i<0 or i>=self.n for i in indices):
            raise ValueError(f'Provide {self.k} indices between 0 and {self.n-1}')
        if len(set(indices))<len(indices):return self.ring.zero()
        inversions=sum(indices[i]>indices[j] for i in range(len(indices)) for j in range(i+1,len(indices)))
        j=self._model['plucker_indices'].index(sorted(indices))
        return (-1)**inversions*self.ring.gen(j)

    def __repr__(self):
        return f'Hypersurface in Gr({self.k},{self.n}), degree {self.d}, {len(self._model["terms"])} terms; use .f for the equation'

    def _latex_(self):
        from .models import latex_equation
        return latex_equation(self._model)

    def _repr_latex_(self):return '$'+self._latex_()+'$'

    def save(self,path):
        """Save the equation; subsequent checks use persistent sidecar checkpoints."""
        path=Path(path).resolve()
        old=self._path
        # Preserve exact bytes: checkpoint identities include the file hash.
        source_matches=(old is not None and old.exists()
                        and json.loads(old.read_text())==self._model)
        if source_matches:
            import shutil
            if old!=path:
                path.parent.mkdir(parents=True,exist_ok=True)
                shutil.copyfile(old,path)
        else:
            write_json(path,self._model)
        if source_matches and old!=path:
            for checkpoint in old.parent.glob(old.name+'.smooth-p*.json'):
                suffix=checkpoint.name[len(old.name):]
                target=path.with_name(path.name+suffix)
                if not target.exists():
                    report=json.loads(checkpoint.read_text())
                    report['model']=str(path)
                    write_json(target,report)
        self._path=path
        return path

    def _model_path(self):
        if self._path is None:
            self._scratch=tempfile.TemporaryDirectory(prefix='grass-sage-')
            self._path=Path(self._scratch.name)/'equation.json'
            write_json(self._path,self._model)
        # Refuse stale external edits: self.f must be the equation actually checked.
        if json.loads(self._path.read_text())!=self._model:
            raise ValueError('The saved equation changed; reload it before checking')
        return self._path

    def check_smoothness(self,*,prime=0,seconds=120,max_charts=1000000,restart=False,verbose=False,
                         backend="persistent",workers=1,algorithm="auto",method="jacobian",
                         sparse=True,symmetry=False,symmetry_permutations=None):
        """Bounded, resumable all-chart test. .is_smooth is True, False, or None.

        No filenames are required. Unsaved objects keep temporary checkpoints
        for their lifetime; call .save(path) to keep them across sessions.
        Reuses one worker by default. Options: workers, algorithm (auto/std/slimgb),
        method (jacobian/contraction for k=3,d=1), symmetry, sparse, backend.
        Successful checkpoints are reused across methods; restart=True recomputes.
        """
        from .smoothness import smooth
        if seconds<=0 or max_charts<0:raise ValueError('seconds must be positive; max_charts nonnegative')
        args=SimpleNamespace(model=str(self._model_path()),prime=int(prime),seconds=float(seconds),
            max_charts=int(max_charts),restart=bool(restart),verbose=verbose,output=None,
            backend=backend,workers=workers,algorithm=algorithm,method=method,dense=not sparse,
            symmetry=symmetry,symmetry_permutations=symmetry_permutations)
        with _messages(verbose):report=smooth(args)
        self.last_smoothness=CheckResult(report,'smoothness')
        return self.last_smoothness

    def check_faithfulness(self,*,max_charts=1000000):
        from .faithfulness import check_faithfulness
        if max_charts<0:raise ValueError('max_charts must be nonnegative')
        self.last_faithfulness=CheckResult(check_faithfulness(self._model,max_charts),'faithfulness')
        return self.last_faithfulness


def load_hypersurface(path):
    """Load an existing CLI or Python equation and reuse its saved checkpoints."""
    from .models import read_model
    return Hypersurface(read_model(path),path)


def load_search(path):
    """Load a saved CLI or Python search as an indexable SearchResults."""
    return SearchResults(json.loads(Path(path).read_text()))
