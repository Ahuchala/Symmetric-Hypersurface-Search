"""Run locally: sage -python output/a7-sparse-form/verify.py.

Checks an exact change of basis from the repository's A7 matrices to the
12-term invariant. No floating point arithmetic or modular inference is used.
The tensor convention is the coefficient representation used by
A7Intersection.M2: the action is exteriorPower(3, r), and similarly for s.
"""
from sage.all import *
from pathlib import Path
import hashlib
import itertools
import json
import re

folder = Path(__file__).resolve().parent
repo = folder.parent.parent
cert = json.loads((folder / 'certificate.json').read_text())
source = (repo / cert['source']).read_bytes()
assert hashlib.sha256(source).hexdigest() == cert['source_sha256']
K = CyclotomicField(7)
zeta7 = K.gen()
b = zeta7 + zeta7**2 + zeta7**4
B = -1-b
matrices = []
for name in ('r', 's'):
    text = re.search(name + r' = matrix\(K, (\{.*?\})\);', source.decode()).group(1)
    matrices.append(matrix(K, sage_eval(text.replace('{', '[').replace('}', ']'), locals={'b': b, 'B': B})))
r, s = matrices
C = matrix(K, [[sage_eval(x, locals={'zeta7': zeta7}) for x in row] for row in cert['change_of_basis']])
assert C.det() != 0
assert r**3 == 1 and s**5 == 1 and (r*s)**7 == 1
triples = list(itertools.combinations(range(10), 3))
terms = [(tuple(t['indices']), K(t['coefficient'])) for t in cert['terms']]

def wedge_image(M):
    return vector(K, [sum(c * M.matrix_from_rows_and_columns(I, J).det() for J, c in terms) for I in triples])

original_tensor = wedge_image(C)
assert original_tensor != 0
assert wedge_image(r*C) == original_tensor
assert wedge_image(s*C) == original_tensor

# Check the displayed cyclic formula, including every exterior-algebra sign.
expected = {}
def add(indices, coefficient):
    inversions = sum(indices[i] > indices[j] for i in range(3) for j in range(i+1, 3))
    key = tuple(sorted(indices))
    expected[key] = expected.get(key, 0) + (-1)**inversions * coefficient
for i in range(3):
    add((0, 1+i, 7+i), 1)
    add((1+i, 4+i, 7+(i+2)%3), 2)
    add((1+i, 1+(i+1)%3, 4+(i+2)%3), -5)
add((1, 2, 3), 9)
add((4, 5, 6), 4)
add((7, 8, 9), 1)
assert dict(terms) == expected

# The same 7-cycle reduces the full invariant-space calculation to 18 columns.
g = C.inverse() * r*s * C
assert g.is_diagonal()
allowed = [J for J in triples if prod(g[j,j] for j in J) == 1]
assert len(allowed) == 18
blocks = []
for M in (C.inverse()*r*C, C.inverse()*s*C):
    blocks.append(matrix(K, [[M.matrix_from_rows_and_columns(I,J).det() - (1 if I == J else 0) for J in allowed] for I in triples]))
assert blocks[0].stack(blocks[1]).right_kernel().dimension() == 1
print('PASS: invertible exact basis change; full A7 invariance; invariant dimension 1; displayed cyclic formula; 12 nonzero integer coefficients.')
print('Basis:', ', '.join(cert['basis']))
print('7-cycle weights:', [next(j for j in range(7) if g[i,i] == zeta7**j) for i in range(10)])
print('Coefficients:', cert['terms'])
