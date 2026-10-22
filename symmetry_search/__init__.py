"""Finite symmetries of Grassmannian hypersurfaces, for interactive Sage.

    from symmetry_search import search, load_hypersurface
    results = search(3, 10, 1, group='A7')
    X = results[0].construct()
    X.f
    X.check_smoothness(prime=113)
"""
from .api import (search, load_hypersurface, load_search, SearchResults,
                  Candidate, Hypersurface, CheckResult)

__all__ = ['search','load_hypersurface','load_search','SearchResults',
           'Candidate','Hypersurface','CheckResult']
