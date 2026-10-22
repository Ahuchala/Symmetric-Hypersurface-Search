"""Necessary conditions for SMOOTH complex DV hyperplane sections.

These are conditional obstructions, not smoothness tests. In particular they
must never prune searches for possibly singular sections. The public switch
smooth_obstructions=False bypasses this module as well as group-order filters.
"""
from sage.all import libgap

SONG_REFERENCE = 'https://arxiv.org/html/2106.13287v2 (Proposition 6.1, Lemma 6.2)'
HM_REFERENCE = 'https://arxiv.org/html/1409.6055 (Theorem 4.9)'

# Hohn-Mason's (order, trace on Leech representation, trace of the square).
# These are the 15 admissible types in Theorem 4.9, NOT traces on H^2 itself.
# The Leech character is H^2 + 1. In the symplectic case H^2 = Q + 3, hence
# the offset from the polarized deformation character Q is FOUR, not three.
DV_SYMPLECTIC_TYPES = (
    (1,24,24), (2,8,24), (3,6,6), (3,-3,-3), (4,4,8),
    (5,4,4), (6,2,6), (6,5,-3), (7,3,3), (8,2,4),
    (9,3,3), (11,2,2), (12,1,5), (14,1,3), (15,1,1),
)


def symplectic_trace_obstruction(tbl, q):
    """Check all classes AND their squares; caller must prove symplecticity."""
    squares = [int(i)-1 for i in libgap.PowerMap(tbl,2)]
    for c, order in enumerate(tbl.OrdersClassRepresentatives()):
        order = int(order)
        allowed = [(t,t2) for n,t,t2 in DV_SYMPLECTIC_TYPES if n == order]
        t, t2 = q[c]+4, q[squares[c]]+4
        if not any(t == a and t2 == b for a,b in allowed):
            return {
                'rule':'dv_symplectic_cohomology_trace',
                'class_index':c+1, 'element_order':order,
                'square_class_index':squares[c]+1,
                'deformation_trace':str(q[c]),
                'deformation_square_trace':str(q[squares[c]]),
                'leech_trace':str(t), 'leech_square_trace':str(t2),
                'allowed_leech_trace_pairs':[list(pair) for pair in allowed],
                'reason':'The perfect group must act symplectically, but the '
                         'predicted cohomological traces are not an admissible '
                         'K3^[2] type.',
                'reference':HM_REFERENCE,
            }
    return None


def smooth_character_obstruction(tbl, chi, sections, psi, k, n, d):
    """Return exact rejection evidence, or None (NOT an existence result).

    Scope: Gr(3,10) hyperplanes, and Gr(7,10) by Grassmannian duality.
    For a smooth section, the projective stabilizer of its equation is finite
    (Song, Prop. 6.1). The derivative of the GL(V) action therefore injects
    End(V) into W tensor psi^{-1}; the scalar direction maps nontrivially to
    the equation. Here W is the 120-dimensional space of equations and psi
    is the character of the selected equation line. Thus

        Q = W tensor psi^{-1} - End(V)

    is an actual 20-dimensional representation. Every irreducible multiplicity
    must be nonnegative. In particular dim(W_psi) >= dim End_G(V), which is
    the centralizer bound. This argument applies to ALL finite groups,
    including cyclic groups, and to both invariant and semi-invariant lines.

    The engine uses S_(1^k)V as its section-character convention. Replacing
    this by the geometric dual convention conjugates BOTH W and psi, hence
    conjugates Q (End(V) is self-dual). Positivity is unchanged, as are the
    rational traces required below. For k=7, det twists cancel in W/psi and
    Grassmannian duality identifies the deformation problem with k=3.

    By Song, Lemma 6.2, Q is the polarized deformation representation. If the
    action is symplectic, contraction with the two-form gives H^2 = Q + 3.
    We apply the trace test ONLY when the group has exactly one linear
    character, equivalently trivial abelianization (a perfect group). This
    includes nonabelian simple groups but does not presume that an arbitrary
    cyclic group is symplectic. Computing the actual symplectic character for
    nonperfect groups is left to a future, separately justified extension.
    """
    if n != 10 or k not in (3,7) or d != 1:
        return None
    q = sections*libgap.ComplexConjugate(psi) - chi*libgap.ComplexConjugate(chi)
    irr = list(tbl.Irr())
    negative = []
    for i, irrep in enumerate(irr):
        multiplicity = int(libgap.ScalarProduct(tbl,q,irrep))
        if multiplicity < 0:
            negative.append({'irreducible_index':i+1,
                             'irreducible_degree':int(irrep[0]),
                             'multiplicity':multiplicity})
    if negative:
        return {
            'rule':'dv_deformation_character_negative_multiplicity',
            'negative_multiplicities':negative,
            'deformation_character_values':[str(v) for v in q],
            'equation_space_dimension':int(libgap.ScalarProduct(tbl,sections,psi)),
            'centralizer_gl_dimension':int(libgap.ScalarProduct(tbl,chi,chi)),
            'reason':'No smooth member: finite projective stabilizer would '
                     'require End(V) to embed in W tensor psi^{-1}, but the '
                     'quotient has a negative irreducible multiplicity.',
            'reference':SONG_REFERENCE,
        }
    if sum(int(irrep[0]) == 1 for irrep in irr) == 1:
        return symplectic_trace_obstruction(tbl,q)
    return None
