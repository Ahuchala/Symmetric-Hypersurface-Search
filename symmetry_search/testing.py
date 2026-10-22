"""Run the full fast correctness suite from Sage or a notebook."""
import unittest


def run_tests(verbosity=1):
    """Return unittest.TestResult; inspect .wasSuccessful(). No benchmarks run."""
    suite=unittest.defaultTestLoader.loadTestsFromNames([
        'symmetry_search.test_workflow',
        'symmetry_search.test_api',
        'symmetry_search.test_regression',
        'symmetry_search.test_optimization',
        'symmetry_search.test_quadratics',
        'symmetry_search.test_smooth_search',
    ])
    return unittest.TextTestRunner(verbosity=verbosity).run(suite)


if __name__=='__main__':
    raise SystemExit(0 if run_tests(verbosity=2).wasSuccessful() else 1)
