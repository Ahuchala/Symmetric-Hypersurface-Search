"""Exact GAP checks of the theorem's scope and search integration."""
import unittest
from unittest.mock import patch
from sage.all import libgap
from symmetry_search import search, load_search
from symmetry_search.search import initialize, smooth_group_obstruction
from symmetry_search.obstructions import smooth_character_obstruction, symplectic_trace_obstruction
import tempfile
from pathlib import Path

class SmoothObstructions(unittest.TestCase):
    @classmethod
    def setUpClass(cls): initialize()

    def test_prime_and_composite_order_exclusions(self):
        self.assertIn(13, smooth_group_obstruction(libgap.CharacterTable('L2(13)'),3,10,1)['forbidden_element_orders'])
        # A9 has only allowed prime divisors, but has order-10 elements.
        self.assertIn(10, smooth_group_obstruction(libgap.CharacterTable('A9'),3,10,1)['forbidden_element_orders'])
        for name in ('A5','A6','A7','L2(11)'):
            self.assertIsNone(smooth_group_obstruction(libgap.CharacterTable(name),3,10,1))

    def test_scope_and_cyclic_exception(self):
        a9=libgap.CharacterTable('A9')
        for k,n,d in ((3,10,2),(3,9,1),(2,10,1)):
            self.assertIsNone(smooth_group_obstruction(a9,k,n,d))
        self.assertIsNotNone(smooth_group_obstruction(a9,7,10,1))
        for p in (2,13,29):
            self.assertIsNone(smooth_group_obstruction(libgap.CyclicGroup(p).CharacterTable(),3,10,1))
        self.assertIsNone(smooth_group_obstruction(libgap.SymmetricGroup(8).CharacterTable(),3,10,1))

    def test_report_optout_and_roundtrip(self):
        with patch('symmetry_search.search.candidates_for_table',return_value=iter(())) as enumerate_reps:
            result=search(3,10,1,group='A9',max_order=1000000)
            enumerate_reps.assert_not_called()
            self.assertTrue(result.complete)
            self.assertFalse(result.errors)
            self.assertEqual(result.report['excluded_groups'][0]['group'],'table:A9')
            self.assertIn('excluded',repr(result))
            with tempfile.TemporaryDirectory() as tmp:
                path=Path(tmp)/'search.json';result.save(path)
                self.assertEqual(load_search(path).report,result.report)
            raw=search(3,10,1,group='A9',max_order=1000000,smooth_obstructions=False)
            enumerate_reps.assert_called_once()
            self.assertEqual(raw.report['excluded_groups'],[])
        with self.assertRaises(ValueError):
            search(3,10,1,smooth_obstructions='yes')

    def test_c2_survives_real_search(self):
        filtered=search(3,10,1,group='small:2:1',reducible=True)
        raw=search(3,10,1,group='small:2:1',reducible=True,smooth_obstructions=False)
        self.assertGreater(len(filtered),0)
        # Cyclic groups are not assumed symplectic, but positivity is necessary
        # for smoothness. Raw mode must retain every excluded family unchanged.
        kept={c.id for c in filtered}
        excluded={x['candidate']['id'] for x in filtered.report['excluded_candidates']}
        self.assertFalse(kept & excluded)
        self.assertEqual(kept | excluded,{c.id for c in raw})
        self.assertTrue(excluded)
        self.assertTrue(all(x['obstruction']['rule']=='dv_deformation_character_negative_multiplicity'
                            for x in filtered.report['excluded_candidates']))
        self.assertEqual(raw.report['excluded_candidates'],[])
        # The invariant 4+6 splitting discussed in chat must remain eligible.
        self.assertTrue(any(c.data['constituents']==[1]*4+[2]*6 and c.multiplicity==64
                            for c in filtered))
        self.assertEqual(filtered.report['excluded_groups'],[])

if __name__=='__main__': unittest.main()
