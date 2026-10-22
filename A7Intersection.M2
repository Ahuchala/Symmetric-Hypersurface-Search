-- A7Intersection.M2
--
-- Computes the A_7-invariant subspace of the 10-dimensional representation
-- on Lambda^3(V) attached to the A_7-symmetric example in Gr(3, 10).
-- The two generators r, s satisfy r^3 = s^5 = id (Atlas presentation),
-- and the invariant subspace is the kernel of (Lambda^3 r - 1, Lambda^3 s - 1).
-- The first column of Wgen gives the coefficient vector fCoeffs of the
-- A_7-invariant linear form f (over QQ(b) with b^2 + b + 2 = 0).
--
-- Run with:  M2 A7Intersection.M2
-- Companion: GriffithsRingA7.M2, A7IrrepDecomposition.M2
-- Thesis reference: Chapter 5, A_7-symmetric example in Gr(3, 10).

needsPackage "Cyclotomic";

K = toField(QQ[t]/(t^2 + t + 2));
b = K_0;
B = -1 - b;

r = matrix(K, {{0,1,0,0,0,0,0,0,0,0},{0,0,0,1,0,0,0,0,0,0},{0,0,0,0,0,1,0,0,0,0},{1,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,1,0},{1,-1/2*B,-b-1/2*B,1/2*B,0,1/2*b-1/4*B,-1/2*B,-1/2*b-1/4*B,-1/2*B,-b-1/2*B},{1,-1/2*B,-1/2*B,-b-1/2*B,1,3/2*b+3/4*B,-1/2*B,1/2*b+3/4*B,-1/2*B,-b-1/2*B},{0,B,0,b,-1,1,-1,-B,-b,b},{-2*b-B,b+3/2*B,b+3/2*B,b-1/2*B,-B,-1/2*b-3/4*B,-b-1/2*B,-1/2*b+1/4*B,-b-1/2*B,1/2*B},{-b+1/2*B,1/2*b-1/4*B,-1/2*b-3/4*B,-3/2*b-1/4*B,-2*b-3/2*B,3/4*b-3/8*B,-3/2*b-5*B/4,3/4*b+13*B/8,-1/2*b-3/4*B,-1/2*b+5/4*B}});

s = matrix(K, {{0,0,1,0,0,0,0,0,0,0},{0,0,0,0,1,0,0,0,0,0},{0,0,0,0,0,0,1,0,0,0},{0,0,0,0,0,0,0,1,0,0},{0,0,0,0,0,0,0,0,0,1},{-b-1/2*B,1/2*b+1/4*B,1/2*b+1/4*B,-1/2*b-3/4*B,-1/2*B,1/4*b-1/8*B,-1/2*b-3/4*B,-3/4*b-1/8*B,1/2*b+1/4*B,-1/2*b-1/4*B},{1,-1/2*B,-b-1/2*B,1/2*B,0,1/2*b-1/4*B,-1/2*B,-1/2*b-1/4*B,-1/2*B,-b-1/2*B},{b+1/2*B,1/2*b+3/4*B,1/2*b+3/4*B,1/2*b-1/4*B,b+1/2*B,-1/4*b+5/8*B,1/2*b-1/4*B,-5/4*b-11/8*B,1/2*b+3/4*B,-1/2*b-3/4*B},{-2,-B,0,-b,0,0,0,0,0,0},{B,2,B,b+2*B,-b,-b-2*B,-b,-2,1,b+2*B}});

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
print Wgen;
fCoeffs = flatten entries Wgen_{0};