# Finite symmetries of Grassmannian hypersurfaces

This collects the character-search approach of `ZetaFunctions/gap_find.gp`,
the reducible-character enumeration in `check1dinv.gap`, the explicit invariant
calculations in `A7Intersection.M2` and `PSL211Intersection.M2`, and smoothness
checking into one local Sage/GAP workflow. Existing scripts are preserved.

Run commands from the repository root with **SageMath**. Sage supplies GAP,
CTblLib, SmallGroups, and Singular; no container or Macaulay2 process is used.
If your usual Sage cache is not writable, prefix commands with
`DOT_SAGE=/tmp/grassmannian-sage`.

Default output is compact. Add `--verbose` to `search`, `construct`, `faithful`,
or `smooth` for detailed diagnostics. Smoothness checks print occasional progress
and one final summary; all chart details remain in the saved evidence file.

## Interactive Sage / Python interface

Open `grassmannian_symmetry.ipynb` using the **SageMath** kernel, or import
from a Sage session whose working directory is the repository root:

```python
from symmetry_search import search, load_hypersurface

results = search(3, 10, 1, group="A7")
results                              # compact, numbered candidate list
X = results[0].construct()
X.f                                  # actual Sage polynomial
X.p(0, 1, 7)                          # a Plücker coordinate
X.invariant_basis                    # Sage polynomials in the same ring
X.generators                         # exact Sage matrices
X.change_of_basis

check = X.check_smoothness(prime=113)
check                                # compact summary
check.is_smooth                       # True, False, or None (inconclusive)
check.report                         # complete evidence when needed
```

No candidate IDs, filenames, or command-line arguments are required. Methods
are quiet by default; pass `verbose=True` to see progress. Use `simple=True`
for simple groups, `simple=False` for non-simple groups, or `simple=None` for
both. `group` accepts one name or a list; omit it to search the chosen bounded
catalogue. Add `reducible=True` to include sums of irreducibles. A search's
`.complete`, `.errors`, and `.report` retain its precise coverage status.

Your already-completed A7 calculation can be loaded without recomputation:

```python
X = load_hypersurface("output/a7-equation.json")
X.f
X.check_smoothness(prime=113)          # reuses all 120 saved chart checks
```

`X.check_faithfulness()` returns a result with `.is_faithful`. Inconclusive
checks return `None`, never a false certificate. Do not use `if check:`;
use `if check.is_smooth is True:` (or `.is_faithful`) explicitly.

For a family, use `results[0].construct(random_section=True, seed=17)` or
`results[0].construct(coefficients=[1,2])`. The invariant basis and polynomial
are representatives in the ambient Plücker polynomial ring, with invariance
verified modulo the Plücker relations. Changing a polynomial expression does
not silently change the verified hypersurface object.

Smoothness checks run bounded workers and keep temporary checkpoints for the
lifetime of an unsaved object. Calling `X.save("output/my-equation.json")`
preserves the equation and its existing checkpoints for a later session.
`results.save(...)` and `load_search(...)` do the same for searches. Direct
Python construction runs in the current Sage process, with the same dimension
and matrix-size limits as the CLI, but no subprocess wall-clock limit; use the
kernel interrupt to stop it. Smoothness retains its per-chart `seconds` limit.

If importing from another working directory, first add the repository folder
to `sys.path`. Ordinary Python without Sage installed cannot perform the exact
algebra; use a Sage-enabled Python interpreter or the SageMath notebook kernel.

## Search, select, construct, check

For the two examples discussed in the paper:

```sh
sage -python grassmannian_symmetry.py search 3 10 1 \
  --group A7 --group 'L2(11)' --output output/search-gr310.json
```

Each result prints a candidate ID. Select a **representation and scalar
character**, not just a group: a group may have several inequivalent actions
and several spaces of invariant hypersurfaces.

```sh
sage -python grassmannian_symmetry.py construct output/search-gr310.json \
  56e5b2e0649b --output output/a7-equation.json

sage -python grassmannian_symmetry.py show output/a7-equation.json
sage -python grassmannian_symmetry.py show output/a7-equation.json --latex

sage -python grassmannian_symmetry.py smooth output/a7-equation.json \
  --prime 113 --seconds 30
```

Use the ID printed by your search if it differs. The saved examples in
`output/symmetry-examples/` include the search and constructed equations.
The A7 basis is named `(z,a0,a1,a2,b0,b1,b2,c0,c1,c2)` in the original sparse
certificate; subsequent optional diagonal rescaling can change coefficients.
All Plücker indices are **zero-based increasing tuples**. Polynomial monomials
in JSON refer to the explicitly stored `plucker_indices` list; they never rely
on Macaulay2's variable ordering.

## Search options and scope

```sh
# Irreducible representations of simple groups in GAP's table library:
sage -python grassmannian_symmetry.py search 2 5 3 --max-order 100000

# Include simple and non-simple groups in that library:
sage -python grassmannian_symmetry.py search 2 5 3 --groups all

# Only non-simple groups:
sage -python grassmannian_symmetry.py search 2 5 3 --groups nonsimple

# Include reducible representations; use a genuinely bounded group catalogue:
sage -python grassmannian_symmetry.py search 1 3 4 \
  --groups all --reducible --catalogue small --max-order 60

# A particular SmallGroups entry:
sage -python grassmannian_symmetry.py search 1 2 3 \
  --groups all --group small:6:1
```

Defaults: simple groups, irreducible representations, order at most 100000,
at most 1000 examined catalogue entries (including entries rejected by filters). `--reducible` includes irreducibles and all sums
of irreducibles of total dimension n, subject to `--max-representations`.
The SmallGroups catalogue is streamed rather than materialized. Its available
orders depend on the local GAP installation. Use a small order bound first.
Library-marked duplicate tables and the `Alt(n)`/`A<n>` spelling are normalized;
results are labelled by their catalogue entries, not claimed to be an exhaustive
isomorphism classification.

This is a **bounded search in the selected catalogue**, not a classification of
all finite groups. It enumerates ordinary linear representations; projective
representations that do not lift to n-dimensional linear representations are
not enumerated. The JSON records coverage, limits, and any errors. A stopped
or partially failed search is not labelled complete.

A search hit proves a **nonzero semi-invariant section exists** and that the
ambient projective action is faithful. It does not yet certify that the group
acts faithfully on a chosen hypersurface, that the equation is reduced, or
that it is smooth. These are separate subsequent checks.

## What is counted

The degree-d coordinate space is the rectangular Schur representation
`S_(d,...,d)(W)`, with k entries, where W is the n-dimensional coordinate
representation. Its character is computed by Newton identities and the
Jacobi–Trudi determinant. This removes the Plücker relations from
`Sym^d(Lambda^k W)`, avoiding false positives that define the zero section.
For d=1 this recovers the original exterior-power character calculation.

All one-dimensional constituents are considered. A section satisfying
`g.f = psi(g) f` defines a G-invariant hypersurface even when psi is nontrivial.
For a nonabelian simple group, psi is necessarily trivial.

The projective kernel is detected exactly by the character condition
`chi(g) conjugate(chi(g)) = n^2`. Checking only `ker(rho)=1` is insufficient:
a nonidentity scalar matrix acts trivially on the Grassmannian.

## Constructing and simplifying equations

Construction uses a basis of standard Plücker monomials (semistandard tableaux
of rectangular shape). For d>1, the generator equations are evaluated on a
Grassmannian big cell. Restriction to this dense cell is injective for sections,
so this calculates invariance **in the quotient coordinate ring**, not in an
unrelated ambient polynomial space. The computed dimension must agree with
the character search.

The program compares the original basis, cyclic eigenbases, and bases adapted
to a cyclic normalizer found by a bounded exact search. It also tries exact
diagonal rescalings to simplify coefficients. The previously certified A7
basis is included as an additional candidate. It chooses by number of nonzero
terms, then coefficient complexity. This is a heuristic, not a proof of global
minimal support. Every accepted equation is rechecked for exact semi-invariance
under every generator. Changes of basis and generator matrices are saved.

The program reuses the repository's A7 and PSL(2,11) matrices when appropriate.
Other selected irreducibles are constructed with GAP's Dixon algorithm;
reducible representations use block sums. Some character-library entries lack
an available concrete group model, or Dixon's construction may fail. In such
a case construction reports the limitation rather than inventing a polynomial.
A compatible explicit GAP group may be supplied with `--gap-group 'EXPRESSION'`.
This argument is trusted local GAP code. AtlasRep may require its data files
for group models not covered by the direct constructors.

Useful construction controls:

- `--basis-tries 0`: keep the original and known certified bases.
- `--seconds 300`: hard overall time limit; `0` disables it.
- `--max-sections 1000`, `--max-entries 2000000`: explicit construction bounds.
- `--max-conductor 420`: limit cyclotomic extensions used for additional bases.
- `--random-section --seed 17`: choose a reproducible small linear combination
  in the original invariant space. Sparse choices can be singular even if a
  smooth member exists.
- `--coefficients 1,2`: choose a specific combination in that same space.
  The full invariant basis is saved with every model.

Changing the seed or combination and rerunning construction lets you test
several members of a higher-dimensional family. A singular chosen member does
not prove that the entire invariant family is singular.

## Faithfulness on the actual hypersurface

Construction tests every nonidentity conjugacy class for a non-fixed-locus
witness on one affine chart by default. When that succeeds, it proves the
geometric action on the **reduced hypersurface** is faithful, even if the
hypersurface is singular. For each class, fixed-point equations are reduced
modulo the squarefree part of the chart equation. A nonzero remainder gives
an exact characteristic-zero witness. Semi-invariance makes it sufficient to
test conjugacy-class representatives.

If one chart is insufficient, the result is `inconclusive`, not `not_faithful`:

```sh
sage -python grassmannian_symmetry.py faithful output/a7-equation.json
```

This checks further charts and saves evidence separately. A nontrivial kernel
is certified only after all charts are checked. Use `--faithfulness-charts 0`
in construction to defer this work. Faithfulness does not assert that the
specified group is the full automorphism group.

## Smoothness, primes, and checkpoints

On each standard affine Grassmannian chart, the code substitutes the minors
of `[I | X]` into f and checks whether `(f, partial f / partial x_ij)` is the
unit ideal. Every one of the `binomial(n,k)` charts must pass for a global
smoothness certificate. This includes all boundary charts.

```sh
# Exact characteristic-zero check (the default):
sage -python grassmannian_symmetry.py smooth output/equation.json

# A potentially faster proof via good reduction:
sage -python grassmannian_symmetry.py smooth output/equation.json --prime 113

# A bounded pilot; explicitly inconclusive until all charts pass:
sage -python grassmannian_symmetry.py smooth output/equation.json \
  --prime 113 --max-charts 2 --seconds 20
```

A smooth nonzero reduction on all charts certifies characteristic-zero
smoothness: the integral equation defines a flat projective Cartier divisor
in the smooth Grassmannian over the local coefficient ring. A **singular
special fibre is inconclusive about characteristic zero**. Denominators,
coefficient-field embeddings, and the possibility of zero reduction are
checked. Non-rational cyclotomic coefficients currently require a splitting
prime in the prime field. A singular characteristic-zero chart certifies
singularity over the algebraic closure, whether or not it has a rational point.

Worker startup and each chart calculation are bounded by `--seconds`; timed-out workers are killed. Timeout/error/partial coverage
is always reported as inconclusive. Repeating the same command resumes the
checkpoint and retries unfinished charts. Checkpoints are tied to a hash of
the exact model; changing equations requires a new output or `--restart`.

## Verification and source files

```sh
sage -python -m unittest symmetry_search.test_workflow -v
```

The tests cover the Schur character and standard-monomial dimensions, removal
of a spurious invariant Plücker relation, projective scalar kernels,
nontrivial scalar characters, reducible enumeration, actual-hypersurface
faithfulness, smooth and singular sections, bad reduction, zero reduction,
and resumable partial coverage. The A7 and PSL(2,11) examples are also exercised
through the command-line workflow and retain their exact certificates.

- `characters.g`, `search.py`: bounded GAP search.
- `representations.py`: concrete group models and matrices.
- `sections.py`: quotient-ring invariants and sparse basis heuristics.
- `construct.py`, `models.py`: construction, exact verification, serialization.
- `faithfulness.py`: the action on the hypersurface itself.
- `smoothness.py`: bounded all-chart checks.
- `data/a7-sparse-certificate.json`: the exact A7 basis from the preceding calculation.

GAP reference: [IrreducibleRepresentationsDixon](https://gap-system.github.io/gap/doc/ref/chap71_mj.html).

## Correctness suite and performance baseline

From a Sage notebook in this repository:

```python
from symmetry_search.testing import run_tests
result = run_tests()
assert result.wasSuccessful()

from symmetry_search.benchmark import run_benchmarks, compare_reports
baseline = run_benchmarks(output="output/my-symmetry-baseline.json")
assert baseline["passed"]
```

The correctness suite covers character search, explicit construction, the public
Python API, exact smoothness and singularity, a singular point visible only on
a boundary chart, nonreduced equations, independent finite-field point checks,
coefficient reduction, worker failures/timeouts, and checkpoint resumption.

A saved reference run is in `symmetry_search/benchmarks/baseline-quick.json`.
After changing an algorithm, run:

```python
trial = run_benchmarks(output="output/symmetry-trial.json")
comparison = compare_reports(trial, baseline)
# Each row has calculation_speedup and wall_speedup; greater than 1 is faster.
```

To compare against the saved reference instead, load it with `json.load`, or
pass `baseline="symmetry_search/benchmarks/baseline-quick.json"` to
`run_benchmarks`. Different Sage/Python/platform versions, fixtures, repeat
counts, or timeout limits are rejected; create your own baseline when these
change. All correctness checks must pass before timings are compared. Code
hashes and Git revision are recorded so modified code can be identified even
before it is committed. Run timings while other expensive computations are
idle, and compare the individual samples as well as medians. Timing ratios
are measurements, not fixed CI pass/fail thresholds.

The default **quick** benchmark has 15 chart cases and three fresh workers per
case: six charts of a smooth Gr(2,4) hyperplane, a boundary singularity, a conic
in characteristic zero and at its bad prime, a nonreduced point, an empty chart,
zero reduction, and three A7 charts at 113. It measures the complete
`check_chart` computation, including chart substitution and ideal membership.
Worker wall time includes Sage startup; `calculation_seconds` excludes startup.
No existing smoothness checkpoint is consulted. Raw results and all timing
samples are retained, including errors and timeouts. A timeout is a failed
benchmark, never a speedup.

```python
full = run_benchmarks(profile="full", output="output/symmetry-full.json")
```

The **full** profile replaces the three A7 charts with all 120 and takes
longer. The quick profile does **not** certify global A7 smoothness; the
previously saved all-chart certificate remains separate. These timings cover
individual chart workers, not parallel scheduling or checkpoint I/O; use the
correctness suite to guard those control paths when optimizing them.

Optional command-line equivalents:

```sh
sage -python -m symmetry_search.testing
sage -python -m symmetry_search.benchmark --output output/symmetry-trial.json
```

## Optimized smoothness checks

Restart an existing Sage notebook kernel after updating the code. The default
now computes only the Plücker minors used by the equation and reuses one Sage
worker across charts. The worker retains the decoded equation and coefficient
field. Existing saved certificates are still reusable.

```python
# Fast general path: a reusable worker and sparse minor construction.
check = X.check_smoothness(prime=113)

# Optional exact chart-orbit reduction and a different Singular algorithm.
check = X.check_smoothness(prime=113, algorithm="slimgb", symmetry=True)

# Two independent Sage workers; useful when individual charts are expensive.
check = X.check_smoothness(prime=113, workers=2)

# For k=3, d=1 only: quadratic contraction equations for the three-form.
check = X.check_smoothness(prime=113, method="contraction", algorithm="slimgb")

# Original computation strategy, retained for comparisons.
check = X.check_smoothness(prime=113, backend="fresh", sparse=False)
```

These calls reuse prior successful chart checks. Pass `restart=True` to
recompute; that replaces the checkpoint for this equation and prime. For timing
comparisons without touching existing evidence, use the benchmark below.
`algorithm` accepts `"auto"` (the original ideal-membership method), `"std"`, or
`"slimgb"`; the algorithm and contraction method remain explicit choices.
`workers` controls separate processes, not concurrent calls into one Singular
instance. Each worker startup and each chart calculation has its own `seconds`
time limit. A timed-out/crashed worker is killed and replaced for subsequent
charts; its unfinished chart is retained as inconclusive and retried on resume.
A batch already running when a singularity is found is allowed to finish.

`method="contraction"` checks all equations
`sigma(u_i,u_j,e_l)=0`, for the three row pairs and all coordinate vectors.
These are quadratic or lower on a chart. They are equivalent to the equation
and its tangent derivatives, including in characteristics 2 and 3. Other
values of `(k,d)` are rejected rather than silently using a different method.

With `symmetry=True`, the program checks monomial saved generators and, for
sparse linear sections with `n <= 12` and at most 100 terms, proposes coordinate
permutations from the support graph. Every accepted proposal must satisfy an
exact ambient-polynomial semi-invariance identity over the checking field,
with all diagonal factors nonzero. This is a sufficient, conservative check;
it can miss symmetries that require Plücker relations or a new basis. It never
assumes the full group permutes the standard charts. The saved A7 equation has
42 verified chart orbits. Reports include the generating permutations, weights,
scalar factors, and the checked representative for each inferred chart.

Explicit coordinate permutations can also be supplied as zero-based image
lists via `symmetry_permutations=[[...], ...]`; an invalid supplied symmetry
raises an error. Finite-field symmetries are used only to cover charts of that
fibre. The resulting all-chart smoothness conclusion uses the same good-
reduction argument as before. `max_charts` limits new representative
computations, so symmetry can certify more charts than that number.

### End-to-end benchmark

The original `run_benchmarks` deliberately starts a fresh process for each
sample and remains useful for chart-algebra comparisons. To measure worker
reuse, parallelism, and chart reduction together with startup and checkpoint
costs, use:

```python
from symmetry_search.pipeline_benchmark import run_pipeline_benchmarks

trial = run_pipeline_benchmarks(
    repeats=3, max_charts=12, output="output/symmetry-pipeline-trial.json"
)
[(v["name"], v["median_wall_seconds"], v.get("speedup"))
 for v in trial["variants"]]

# Full A7 verification, including the symmetry-reduced variant:
full = run_pipeline_benchmarks(
    repeats=1, max_charts=None, output="output/symmetry-pipeline-full.json"
)
```

The default fixture is the known smooth A7 equation at 113. Each sample starts
with a fresh checkpoint. The variants are the dense/fresh reference,
persistent/sparse, two-worker parallel, `slimgb`, quadratic contractions, and
(for full runs) verified symmetry reduction. Three-form contractions are
omitted for other `(k,d)`. The runner accepts a different `model` path and
`prime`, intended for other known smooth fixtures. Partial runs must have the
specified chart coverage and an incomplete result; full runs must certify
smoothness. Timing ratios are reported only if every sample passes. Reports
retain complete per-chart evidence, execution settings, code hashes, and
software versions. The saved full comparison is
`benchmarks/pipeline-a7-full.json`; its timings are one sample per variant,
not repeated statistical estimates.

## Quadratic construction for general Gr(k,n)

Degree-two candidates now use an exact Reynolds-form constructor automatically:

```python
results = search(3, 10, 2, group="A6")
X = results[0].construct(verbose=True)
len(X.invariant_basis)                 # 25
X.data["quadratic_construction"]       # dimensions and exact verification evidence
```

The method works for arbitrary `(k,n)` in characteristic zero, including
reducible representations and nontrivial scalar characters. Resource bounds
still apply; this is not a claim that every large group or exterior power is
cheap. The Plücker representation has dimension `binomial(n,k)`.

For a scalar character `psi`, the constructor averages rank-one symmetric
forms using `sum_g psi(g)^(-1) (g*v)(g*v)^t`. It verifies each accepted form's
semi-invariance exactly under the generators. An exact group enumeration
using the original n-dimensional matrices supplies the words for averaging;
its size must equal the certified group order. For cyclotomic representations,
the actual coefficient subfield is found by exact rational linear algebra
and closure under multiplication. Orbit propagation tracks integer numerators and exact common denominators;
most Gram-matrix work also uses the native integer/rational matrix backend.
This avoids generic number-field arithmetic without using approximations.

The resulting quadrics are straightened to standard Plücker monomials using
integer polynomial identities in a generic k-by-n matrix. Thus Plücker
relations do not count as equations. Exact independent evaluations on
explicit decomposable vectors prove independence in the Grassmannian
coordinate ring. Construction stops only when that rank equals the
character-search multiplicity, proving completeness. Random seeds assist
discovery; every accepted result is certified exactly. A bounded unsuccessful
search raises `LimitExceeded`, never returns a partial space as complete.

RREF makes the final invariant basis canonical in the original coordinates
before any optional diagonal rescaling. `coefficients=[...]` refers to this
basis; `random_section=True` chooses a reproducible small combination. Sparse
choices can still be singular or reducible. Smoothness and faithfulness remain
separate checks. The quadratic path uses RREF and diagonal simplification;
`basis_tries` and `max_conductor` control only the older expanded constructor.

```python
# Retain the old algorithm for small-case comparisons:
Y = results[0].construct(quadratic_method="expanded", max_sections=5000)

# Bounds for the new method:
X = results[0].construct(max_entries=2_000_000, quadratic_trials=120)
```

`max_sections=None` is now the default: expanded construction retains its
1000-section cap, while the quadratic method checks its orbit, generator,
straightening-cache, and output-basis entry budgets through `max_entries`.
An explicitly supplied `max_sections` is always respected. `quadratic_trials=0`
selects a dimension-dependent trial bound; change the seed or increase the
bound if exact quotient independence is not reached. Entry bounds are not
byte-accurate RAM estimates. Direct notebook construction is interruptible
but has no hard wall-clock limit; the CLI retains `--seconds`.

The tests compare with the expanded algorithm, check all quadratic
straightening identities in several small Grassmannians, test scalar
characters and cyclotomic fields, and verify budget failures. The saved
`output/symmetry-examples/a6-gr310-quadric.json` contains the full 25-dimensional
basis, chosen equation, matrices, and verification evidence.
The corresponding seven-dimensional A7 example is saved as
`output/symmetry-examples/a7-gr310-quadric.json`. Neither example has a
smoothness certificate. Single-run timings and source/model hashes are in
`benchmarks/quadratic-pilots.json` (A6 about 20 seconds including search and
faithfulness; A7 about 127 seconds including process startup, with faithfulness
skipped). These runs overlapped other tests and are not isolated benchmarks.

## Immediate singularity detection

Before launching a Gröbner worker, smoothness now checks coordinate planes.
For a nonzero standard-monomial equation, at the plane with `p_I=1`, its value
is the coefficient of `p_I^d`. Tangent derivatives are, up to sign, the
coefficients of `p_I^(d-1)*p_J` for `|I intersection J|=k-1`. Their simultaneous
vanishing provides an exact singular point. Zero reductions and nonstandard
representatives are conservatively left to the existing chart engine.

This detects the saved A5 cubic on Gr(2,5) as **singular in characteristic
zero** without a Gröbner basis. A modular singular witness still establishes
only singularity of that special fibre. The precheck is skipped when
`max_charts=0`, which remains a request to inspect existing checkpoints only.
For other difficult equations, `seconds=10, max_charts=1, verbose=True` gives
a short visible pilot; the default timeout is per chart, not per entire run.

## Continue until a smooth, faithful example is certified

```python
results = search(3, 10, 2)
X = results.find_smooth(
    sections_per_candidate=5,  # sparse equation, then up to four combinations
    prime=113, seconds=10,    # timeout per smoothness chart
    algorithm="slimgb", symmetry=True,
    seed=0, verbose=True,
)
if X is not None:
    X.f
    X.last_smoothness
    X.last_faithfulness

results.last_smooth_search  # outcomes and evidence for all attempts
```

Candidates are tried in result order. Each invariant space is constructed
once, then reused for reproducible random combinations. A one-dimensional
space has only one equation up to scale, so it is tried only once. Duplicate
scalar multiples are skipped. `sections_per_candidate` bounds the total
number of equations per candidate; `max_candidates` optionally limits how
many candidate rows are visited. Distinct representations and scalar
characters of the same group remain separate candidates.

Success requires **both** an exact faithful-action certificate on the
hypersurface and a characteristic-zero smoothness certificate. The default
`prime=113` uses smooth good reduction to certify characteristic-zero
smoothness; `prime=0` checks directly. A singular special fibre, a timeout,
or an inconclusive faithfulness check never counts as success. Such outcomes
are recorded and the search continues. Construction limits/errors are also
recorded; programming errors and user interrupts are not silently swallowed.

`None` means no certificate was obtained within the chosen bounds, **not**
that the groups admit no smooth invariant hypersurface. A bad reduction
prime can hide a smooth characteristic-zero example; inspect the report and
retry with another prime or `prime=0`. There is no overall time limit:
`seconds` is per smoothness chart, while construction and faithfulness run
in the calling Sage process and can be interrupted. `max_charts` and
`faithfulness_charts` bound chart coverage separately.

Pass constructor limits through `construction_options`, for example
`construction_options={"max_entries": 4000000}`. The method controls the seed,
section choices, construction verbosity, and separate faithfulness check.
Restart the Sage kernel after updating the package.

## Smooth searches and searches allowing singular sections

The switch `smooth_obstructions` controls **all** pruning based on smoothness.
For future work on singular sections, use:

```python
# Includes both smooth and singular possibilities; not a singularity certificate.
results = search(3, 10, 1, reducible=True, smooth_obstructions=False)
```

The CLI equivalent is `--no-smooth-obstructions`. This recovers the raw
representation/section-character candidates, including their original IDs.
It still enforces the ordinary search requirements (e.g. a faithful projective
ambient representation and a nonzero semi-invariant equation space).
Construction, faithfulness checks, and saving/loading remain available in
this mode. `find_smooth` is specifically for smooth examples.

The default `smooth_obstructions=True` applies these **necessary**, not
sufficient, conditions to degree-one sections on `Gr(3,10)` and `Gr(7,10)`:

1. A nonabelian simple group cannot have element orders outside
   `{1,2,3,4,5,6,7,8,9,11,12,14,15}`. Rejections are recorded in
   `results.report["excluded_groups"]`.
2. For every finite group, put `Q = W tensor psi^{-1} - End(V)`, where `W`
   is the section representation and `psi` the equation-line character.
   Every irreducible multiplicity of this virtual character must be
   nonnegative. In particular `dim W_psi >= dim End_G(V)`.
3. For perfect groups (those with no nontrivial linear characters), the
   induced action on the Debarre-Voisin fourfold must be symplectic.
   The predicted Leech character is `Q + 4`; its traces on each element
   and its square must match Hohn-Mason's admissible types.

Candidate-level exclusions are recorded in
`results.report["excluded_candidates"]`. Each entry contains the original
candidate and exact obstruction evidence. Summaries show exclusion counts;
these are proven exclusions, not computational errors. Saved reports preserve
this evidence. Already-saved searches are not retroactively filtered.

The proofs, dual/section-character conventions, and hypotheses are commented
in `search.py` and `obstructions.py`. Sources:

- [Song, Proposition 6.1 and Lemma 6.2](https://arxiv.org/html/2106.13287v2):
  finite stabilizers and polarized deformations.
- [Hohn-Mason, Theorem 4.9](https://arxiv.org/html/1409.6055): symplectic types.
- [Benedetti-Manivel, Theorem 1.1](https://arxiv.org/abs/2110.00255): extension
  of smooth hyperplane automorphisms to the Grassmannian.

Cyclic groups are **not presumed symplectic**. They are subject to test 2 in
smooth mode, but not test 3. No order-22 cutoff is imposed. These tests do not
apply to other Grassmannians or higher degrees. They do not replace checking
a selected equation for smoothness and hypersurface faithfulness, nor do they
bound group order or fill gaps in the catalogue. Completeness remains relative
to the stated catalogue, ordinary linear representations, and configured bounds.
