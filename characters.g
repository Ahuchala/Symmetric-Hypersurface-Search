# Character calculations for sections of O(d) on Gr(k,n).
# Uses S_(d^k), not Sym^d(Lambda^k): the latter includes Plucker relations.
GrassSectionCharacter := function(tbl, chi, k, d)
    local m, c, j, q, pm, h, vals, mat, maxm;
    maxm := d+k-1;
    pm := List([1..maxm], j -> PowerMap(tbl,j));
    vals := [];
    for c in [1..Length(chi)] do
        h := [1];
        for m in [1..maxm] do
            Add(h, Sum([1..m], j -> h[m-j+1]*chi[pm[j][c]])/m);
        od;
        mat := [];
        for m in [1..k] do
            Add(mat, []);
            for j in [1..k] do
                q := d-m+j;
                if q < 0 then Add(mat[m],0); else Add(mat[m],h[q+1]); fi;
            od;
        od;
        Add(vals,DeterminantMat(mat));
    od;
    return Character(tbl,vals);
end;

GrassProjectiveKernelSize := function(tbl,chi)
    local n, sizes;
    n := chi[1];
    sizes := SizesConjugacyClasses(tbl);
    return Sum(Filtered([1..Length(chi)],i -> chi[i]*ComplexConjugate(chi[i])=n*n),i -> sizes[i]);
end;

GrassCharacterValue := function(G,chi,g)
    return chi[PositionProperty(ConjugacyClasses(G),C -> g in C)];
end;
