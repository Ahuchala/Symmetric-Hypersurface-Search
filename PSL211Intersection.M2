-- PSL211Intersection.M2
--
-- Computes the PSL(2, 11)-invariant subspace of Lambda^3(V_10) for the
-- 10-dimensional rational representation of PSL(2, 11) on V_10. The two
-- generators r, s satisfy the standard PSL(2, 11) presentation; the
-- kernel of (Lambda^3 r - 1, Lambda^3 s - 1) gives the invariant subspace
-- and fCoeffs (= -1 * first column) gives the coefficient vector of the
-- PSL(2, 11)-invariant linear form f used by GriffithsRingPSL211.M2.
--
-- Run with:  M2 PSL211Intersection.M2
-- Companion: GriffithsRingPSL211.M2, PSL211IrrepDecomposition.M2
-- Thesis reference: Chapter 5, PSL(2, 11)-symmetric example in Gr(3, 10).

K = QQ;


r = matrix(K, {
    {-1, 2, 0, 1, 0,-1, 2,-1, 1,-1},
    {-1, 1, 1, 0, 0,-1, 1,-1, 0, 0},
    { 0, 1, 0, 0,-1,-1, 0, 0, 0,-1},
    { 0, 0, 1, 0, 1,-1, 0,-1, 1, 1},
    {-1, 0, 1, 0, 0,-1, 0, 0, 0, 1},
    { 0, 1, 0, 0,-1,-1, 1,-1, 0,-1},
    { 0, 1,-1, 1,-1, 0, 1, 0, 0,-1},
    { 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
    { 1,-1,-1, 0, 0, 1,-1, 0, 0,-1},
    { 0, 0, 0, 0, 0, 1, 0, 0, 0, 0}});

s = matrix(K, {
    {-1, 2, 1, 0, 0,-2, 2,-1, 1, 0},
    { 0, 1, 0, 0,-1,-1, 1,-1, 0,-1},
    { 0, 1, 0, 0, 0, 0, 0, 0, 0, 0},
    { 1,-1,-1,-1, 0, 1,-2, 1, 0, 0},
    { 0, 0, 0, 0, 0,-1, 0,-1, 0, 0},
    { 0, 0, 0, 0, 0, 1, 0, 0, 0, 0},
    {-1, 1, 1, 1, 1,-1, 2,-1, 1, 1},
    {-1, 2, 0, 1, 0,-2, 2,-1, 1, 0},
    {-1, 1, 0, 0,-1, 0, 0, 1, 0, 0},
    { 0, 0, 0, 0, 0, 0, 1, 0, 0, 0}});
n = 10;
triples = flatten flatten for i from 0 to n-1 list
            for j from i+1 to n-1 list
              for l from j+1 to n-1 list {{i,j,l}};

wr = exteriorPower(3,r);
ws = exteriorPower(3,s);
-- I  = id_(K^m);

Aeq  = (wr - 1) || (ws - 1);   -- fixed: || not **
Wgen = mingens ker Aeq;

print("Invariant subspace dimension: " | toString numgens source Wgen);
-- print Wgen;
fCoeffs = -1 * flatten entries Wgen_{0}