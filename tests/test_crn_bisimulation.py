#!/usr/bin/env python
#
#  tests/test_crn_bisimulation.py
#  Original source from the Nuskell compiler project
#
import logging
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

import unittest

from crnverifier.utils import parse_crn
from crnverifier.crn_bisimulation import (crn_bisimulation_test, 
                                          modular_crn_bisimulation_test,
                                          crn_bisimulations)

SKIP_SLOW = True

class HelperTests(unittest.TestCase):
    def test_modularity_condition(self):
        from collections import Counter
        from crnverifier.crn_bisimulation import moduleCond
        module = """ B <=> e101
                     e104 -> e103 + e105
                     e101 -> e102 + e103 + e104
                 """
        module = [[Counter({'B': 1}), Counter({'e101': 1})], 
                  [Counter({'e101': 1}), Counter({'B': 1})], 
                  [Counter({'e104': 1}), Counter({'e103': 1, 'e105': 1})],
                  [Counter({'e101': 1}), Counter({'e102': 1, 'e103': 1, 'e104': 1})]]
        fsc = {'A', 'B', 'C'}
        isc = {'A', 'B', 'C'}

        minter1 = {'A': Counter({'A': 1}), 
                   'B': Counter({'B': 1}), 
                   'C': Counter({'C': 1}), 
                   'e101': Counter({'B': 1}),
                   'e102': Counter(), 
                   'e103': Counter(), 
                   'e104': Counter({'B': 1}),
                   'e105': Counter()} 
        assert moduleCond(module, fsc, isc, minter1) is False

        minter2 = {'A': Counter({'A': 1}), 
                   'B': Counter({'B': 1}), 
                   'C': Counter({'C': 1}), 
                   'e101': Counter({'B': 1}), 
                   'e102': Counter(),
                   'e103': Counter(), 
                   'e104': Counter(), 
                   'e105': Counter()}
        self.assertNotEqual(minter1, minter2)
        assert moduleCond(module, fsc, isc, minter2) is True

        fsc = {'B'}
        isc = {'B'}
        minter3 = {'B': Counter({'B': 1}), 
                   'e101': Counter({'B': 1}), 
                   'e102': Counter(),
                   'e103': Counter(), 
                   'e104': Counter({'B': 1}), 
                   'e105': Counter()}
        self.assertNotEqual(minter1, minter3)
        assert moduleCond(module, fsc, isc, minter3) is False

        minter4 = {'B': Counter({'B': 1}), 
                   'e101': Counter({'B': 1}), 
                   'e102': Counter(),
                   'e103': Counter(), 
                   'e104': Counter(), 
                   'e105': Counter()}
        self.assertNotEqual(minter3, minter4)
        assert moduleCond(module, fsc, isc, minter4) is True

    def test_modularity_condition_02(self):
        from collections import Counter
        from crnverifier.crn_bisimulation import moduleCond
        module = [[Counter({'A': 1}), Counter({'e106': 1})],
                  [Counter({'e106': 1}), Counter({'A': 1})],
                  [Counter({'A': 1}), Counter({'e100': 1})], 
                  [Counter({'e100': 1}), Counter({'A': 1})],
                  [Counter({'e100': 1}), Counter({'e106': 1, 'A': 1})],
                  [Counter({'e106': 1, 'A': 1}), Counter({'e100': 1})]]
        fsc, isc = {'A', 'B', 'C'}, set()
        inter = {'A': Counter({'A': 1}), 
                 'e100': Counter({'A': 2}),
                 'e106': Counter({'A': 1})}
        assert moduleCond(module, fsc, isc, inter) is False
        fsc, isc = {'A'}, {'A'}
        inter  = {'A': Counter({'A': 1}), 
                  'e100': Counter({'A': 2}), 
                  'e106': Counter({'A': 1})}
        assert moduleCond(module, fsc, isc, inter) is True

        #imodule = "a + b -> c + d"
        #fsm = set('ABCD')
        #ism = set('abcd')
        #minter = {'a': 'A',
        #          'b': 'B',
        #          'c': 'C',
        #          'd': 'D'}
        #x = modularity_condition(imodule, fsm, ism, minter)
        #print(x)

class BisimulationTests(unittest.TestCase):
    def test_QingDong_thesis(self):
        # An example where the choice of the permissive checker matters ...
        fcrn, fs = parse_crn('tests/crns/crn6.crn', is_file = True)
        icrn, _ = parse_crn('tests/crns/icrns/crn6_qingdong_thesis.crn', is_file = True)

        inter_01 = {'i778': ['Y'],
                    'i575': ['X'],
                    'i599': ['C'],
                    'i2232': ['A'],
                    'i73': ['B']}

        inter_02 = {'i842': ['Y', 'X', 'A'],
                    'i394': ['X', 'Y', 'X'],
                    'i119': ['X', 'B', 'A'],
                    'i2300': ['A', 'C'],
                    'i778': ['Y'],
                    'i575': ['X'],
                    'i599': ['C'],
                    'i2232': ['A'],
                    'i73': ['B']}

        v, _ = crn_bisimulation_test(fcrn, icrn, fs, 
                                     interpretation = inter_01,
                                     permissive = 'graphsearch')
        self.assertTrue(v)

        v, _ = crn_bisimulation_test(fcrn, icrn, fs, 
                                     interpretation = inter_01, 
                                     permissive = 'reactionsearch')
        self.assertTrue(v)

        if not SKIP_SLOW:
            # TODO: How long approximately?
            v, _ = crn_bisimulation_test(fcrn, icrn, fs, 
                                         interpretation = inter_01, 
                                         permissive = 'loopsearch')
            self.assertTrue(v)

        v, _ = crn_bisimulation_test(fcrn, icrn, fs,
                                     interpretation = inter_02, 
                                     permissive = 'graphsearch')
        self.assertTrue(v)

        v, _ = crn_bisimulation_test(fcrn, icrn, fs, 
                                     interpretation = inter_02, 
                                     permissive = 'reactionsearch')
        self.assertTrue(v)

        v, _ = crn_bisimulation_test(fcrn, icrn, fs,
                                     interpretation = inter_02, 
                                     permissive = 'loopsearch')
        self.assertTrue(v)

        if not SKIP_SLOW:
            # These tests complete in less than 10 minutes
            v, _ = crn_bisimulation_test(fcrn, icrn, fs, permissive = 'graphsearch')
            self.assertTrue(v)

            v, _ = crn_bisimulation_test(fcrn, icrn, fs, permissive = 'reactionsearch')
            self.assertTrue(v)

            # These might not even finish in an overnight run ... who knows?
            #v, _ = crn_bisimulation_test(fcrn, icrn, fs, permissive = 'loopsearch')
            #self.assertTrue(v)

    @unittest.skipIf(SKIP_SLOW, "skipping slow tests")
    def test_qian_roessler_bisimulation(self):
        # TODO: How long approximately?
        (fcrn, fs) = parse_crn('tests/crns/roessler_01.crn', is_file = True)
        (icrn, _) = parse_crn('tests/crns/icrns/roessler_qian2011.crn', is_file = True)

        partial = {sp: [sp] for sp in fs}
        backup = {sp: [sp] for sp in fs}
        v, i = crn_bisimulation_test(fcrn, icrn, fs, partial)
        self.assertTrue(v)
        self.assertDictEqual(partial, backup)

    def test_example_01(self):
        # A sample test to aggree on a new interface for bisimulation.
        fcrn = "A->B"
        ecrn = "A<=>i19; i19<=>i39+X; i39->i71+i72"

        fcrn, fs = parse_crn(fcrn)
        ecrn, _ = parse_crn(ecrn)
        partial = {sp: [sp] for sp in fs}

        v, i = crn_bisimulation_test(fcrn, ecrn, fs, 
                                     interpretation = partial,
                                     permissive = 'graphsearch')
        self.assertTrue(v)

        v, i = crn_bisimulation_test(fcrn, ecrn, fs, 
                                     interpretation = partial,
                                     permissive = 'loopsearch')
        self.assertTrue(v)

        v, i = crn_bisimulation_test(fcrn, ecrn, fs, 
                                     interpretation = partial,
                                     permissive = 'reactionsearch',
                                     permissive_depth = 8)
        self.assertTrue(v)

        # A function that does not say so, should not modify its arguments.
        self.assertDictEqual(partial, {sp: [sp] for sp in fs})

    def test_example_02(self):
        fcrn = """A + B -> C + D
                  A + C -> B + D"""
        icrn = """x1 -> x2
                  x3 + x4 <=> x5
                  x2 -> x6 + x8
                  x5 -> x7
                  x3 <=> x6
                  x9 <=> x10
                  x10 + x4 <=> x1 
                  x7 -> x9 + x8"""

        # First correct interpretation
        inter1 = {'x1': ['A', 'B'], 
                  'x2': ['C', 'D'],
                  'x3': ['C'],
                  'x4': ['A'],
                  'x5': ['A', 'C'],
                  'x6': ['C'],
                  'x7': ['B', 'D'],
                  'x8': ['D'],
                  'x9': ['B'],
                  'x10': ['B']}
        pinter1 = {'x7': ['B', 'D']}

        # Second correct interpretation
        inter2 = {'x1': ['A', 'C'], 
                  'x2': ['B', 'D'],
                  'x3': ['B'],
                  'x4': ['A'],
                  'x5': ['A', 'B'],
                  'x6': ['B'],
                  'x7': ['C', 'D'],
                  'x8': ['D'],
                  'x9': ['C'],
                  'x10': ['C']}
        pinter2 = {'x7': ['C', 'D']}

        # CRN preprocessing
        fcrn, fs = parse_crn(fcrn)
        icrn, _ = parse_crn(icrn)

        # Using partial inter1
        v, i1 = crn_bisimulation_test(fcrn, icrn, fs, 
                                      interpretation = pinter1,
                                      permissive = 'graphsearch')
        self.assertTrue(v)
        self.assertDictEqual(inter1, i1)

        v, i1 = crn_bisimulation_test(fcrn, icrn, fs, 
                                      interpretation = pinter1,
                                      permissive = 'loopsearch')
        self.assertTrue(v)
        self.assertDictEqual(inter1, i1)

        v, i1 = crn_bisimulation_test(fcrn, icrn, fs, 
                                      interpretation = pinter1,
                                      permissive = 'reactionsearch')
        self.assertTrue(v)
        self.assertDictEqual(inter1, i1)

        # Using inter1
        v, i1 = crn_bisimulation_test(fcrn, icrn, fs, 
                                      interpretation = inter1,
                                      permissive = 'graphsearch')
        self.assertTrue(v)
        self.assertDictEqual(inter1, i1)

        v, i1 = crn_bisimulation_test(fcrn, icrn, fs, 
                                      interpretation = inter1,
                                      permissive = 'loopsearch')
        self.assertTrue(v)
        self.assertDictEqual(inter1, i1)

        v, i1 = crn_bisimulation_test(fcrn, icrn, fs, 
                                      interpretation = inter1,
                                      permissive = 'reactionsearch')
        self.assertTrue(v)
        self.assertDictEqual(inter1, i1)

        # Using partial inter2
        v, i2 = crn_bisimulation_test(fcrn, icrn, fs, 
                                      interpretation = pinter2,
                                      permissive = 'graphsearch')
        self.assertTrue(v)
        self.assertDictEqual(inter2, i2)

        v, i2 = crn_bisimulation_test(fcrn, icrn, fs, 
                                      interpretation = pinter2,
                                      permissive = 'loopsearch')
        self.assertTrue(v)
        self.assertDictEqual(inter2, i2)

        v, i2 = crn_bisimulation_test(fcrn, icrn, fs, 
                                      interpretation = pinter2,
                                      permissive = 'reactionsearch')
        self.assertTrue(v)
        self.assertDictEqual(inter2, i2)

    def test_example_02_false(self):
        fcrn = """A + B -> C + D
                  A + C -> B + D"""
        icrn = """x1 -> x2
                  x3 + x4 <=> x5
                  x2 -> x6 + x8
                  x5 -> x7
                  x3 <=> x6
                  x9 <=> x10
                  x10 + x4 <=> x1 
                  x7 -> x9 + x8"""

        fcrn, fs = parse_crn(fcrn)
        icrn, _ = parse_crn(icrn)

        v, _ = crn_bisimulation_test(fcrn, icrn, fs)
        self.assertTrue(v)

        # Test wrong partial interpretation
        partial = dict()
        partial['x2'] = ['B', 'D']
        v, _ = crn_bisimulation_test(fcrn, icrn, fs, 
                                     interpretation = partial,
                                     permissive = 'loopsearch')
        self.assertTrue(v)

        partial['x3'] = ['C']
        v, _ = crn_bisimulation_test(fcrn, icrn, fs, 
                                     interpretation = partial,
                                     permissive = 'loopsearch')
        self.assertFalse(v)

    def test_example_04(self):
        # Two valid interpretations
        fcrn = "B + B -> B"
        icrn = "B <=> x1; B + x1 -> x2 + x3; x2 -> B + x4"

        fcrn, fs = parse_crn(fcrn)
        icrn, _ = parse_crn(icrn)

        ifull1 = {'B': ['B'],
                  'x1': ['B'],
                  'x2': ['B', 'B'],
                  'x3': [],
                  'x4': []}

        ipart1 = {'B': ['B'],
                  'x2': ['B', 'B']}

        v, i1 = crn_bisimulation_test(fcrn, icrn, fs, 
                                      interpretation = ipart1)
        self.assertTrue(v)
        self.assertDictEqual(i1, ifull1)

        ifull2 = {'B': ['B'],
                  'x1': ['B'],
                  'x2': ['B'],
                  'x3': [],
                  'x4': []}

        ipart2 = {'B': ['B'],
                  'x2': ['B']}

        v, i2 = crn_bisimulation_test(fcrn, icrn, fs, 
                                      interpretation = ipart2)
        self.assertTrue(v)
        self.assertDictEqual(i2, ifull2)

    def test_example_05(self):
        # Issue fixed: Naming species in certain ways broke bisimulation
        fcrn = "A+C->A+B"
        icrn = "A <=> x1 + e45; C + x1 <=> x3 + x4; x3 -> A + B + x5"

        fcrn, fs = parse_crn(fcrn)
        icrn, _ = parse_crn(icrn)

        inter = {'A': ['A'],
                 'B': ['B'],
                 'C': ['C']}

        v, i1 = crn_bisimulation_test(fcrn, icrn, fs, interpretation = inter)
        self.assertTrue(v)

    def test_garbage_collection(self):
        # Garbage collection schemes do not produce a correct CRN bisimulation ...
        fcrn = "A + B <=> X + Y"
        icrn = """ 
                A <=> i22
                i59 <=> i139
                i45 -> i351 + i352
                i22 + B <=> i45 + i44
                i44 <=> i60 + i59
                i60 -> i104 + i105
                i139 <=> i227 + X
                i227 <=> i269 + Y
                i269 -> i338 + i339
               """

        fcrn, fs = parse_crn(fcrn)
        icrn, _ = parse_crn(icrn)

        inter = {'A': ['A'],
                 'B': ['B'],
                 'X': ['X'],
                 'Y': ['Y'],
                 'i22': ['A'],
                 'i44': ['A', 'B'],
                 'i59': ['A', 'B'],
                 'i139': ['A', 'B'],
                 'i227': ['Y']}

        v, i1 = crn_bisimulation_test(fcrn, icrn, fs, interpretation=inter)
        assert v is False

    def test_order_formals_bug(self):
        fcrn = "A + B -> C"
        icrn = "B_1_ + i7 -> i684 + i17; A <=> i7; i17 -> C_1_ + i29"

        fcrn, fs = parse_crn(fcrn)
        icrn, _ = parse_crn(icrn)

        print(fs)
        i01 = {'A': ['A'], 'i7': ['A'], 'B_1_': ['B'], 'i684': [],    'i17': ['C'],      'i29': [],    'C_1_': ['C']   }
        i02 = {'A': ['A'], 'i7': ['A'], 'B_1_': ['B'], 'i684': [],    'i17': ['C'],      'i29': ['C'], 'C_1_': []      }
        i03 = {'A': ['A'], 'i7': ['A'], 'B_1_': ['B'], 'i684': [],    'i17': ['A', 'B'], 'i29': [],    'C_1_': ['C']   }
        i04 = {'A': ['A'], 'i7': ['A'], 'B_1_': ['B'], 'i684': [],    'i17': ['A', 'B'], 'i29': ['C'], 'C_1_': []      }
        i05 = {'A': ['A'], 'i7': ['A'], 'B_1_': ['B'], 'i684': ['C'], 'i17': [],         'i29': [],    'C_1_': []      }
        i06 = {'A': ['B'], 'i7': ['B'], 'B_1_': ['A'], 'i684': [],    'i17': ['C'],      'i29': [],    'C_1_': ['C']   }
        i07 = {'A': ['B'], 'i7': ['B'], 'B_1_': ['A'], 'i684': [],    'i17': ['C'],      'i29': ['C'], 'C_1_': []      }
        i08 = {'A': ['B'], 'i7': ['B'], 'B_1_': ['A'], 'i684': [],    'i17': ['A', 'B'], 'i29': [],    'C_1_': ['C']   }
        i09 = {'A': ['B'], 'i7': ['B'], 'B_1_': ['A'], 'i684': [],    'i17': ['A', 'B'], 'i29': ['C'], 'C_1_': []      }
        i10 = {'A': ['B'], 'i7': ['B'], 'B_1_': ['A'], 'i684': ['C'], 'i17': [],         'i29': [],    'C_1_': []      }


        inters = list(crn_bisimulations(fcrn, icrn))[1:]
        for r in inters:
            for e, i in enumerate([i01, i02, i03, i04, i05, i06, i07, i08, i09, i10], 1):
                if r == i:
                    print(e, r)
                    break

        assert all(i in inters for i in [i01, i02, i03, i04, i05, i06, i07, i08, i09, i10])


class ModularBisimulationTests(unittest.TestCase):
    def test_qian_roessler_modular_bisimulation(self):
        (fcrns, fs) = parse_crn('tests/crns/roessler_01.crn', is_file = True, modular = True)
        icrns, _ = parse_crn('tests/crns/icrns/roessler_qian2011_modular.crn', is_file = True, modular = True)
        partial = {sp: [sp] for sp in fs}
        backup = {sp: [sp] for sp in fs}
        v, i = modular_crn_bisimulation_test(fcrns, icrns, fs, partial)
        self.assertTrue(v)
        self.assertDictEqual(partial, backup)
        v, i = modular_crn_bisimulation_test(fcrns, icrns, fs)

    def test_qian_roessler_bisimulation_with_modular_interpretation(self):
        fcrns, fs = parse_crn('tests/crns/roessler_01.crn', is_file = True, modular = True)
        icrns, _ = parse_crn('tests/crns/icrns/roessler_qian2011_modular.crn', is_file = True, modular = True)
        partial = {sp: [sp] for sp in fs}
        v, i = modular_crn_bisimulation_test(fcrns, icrns, fs, partial)
        self.assertTrue(v)

        fcrn, _ = parse_crn('tests/crns/roessler_01.crn', is_file = True)
        icrn, _ = parse_crn('tests/crns/icrns/roessler_qian2011.crn', is_file = True)
        v, i = crn_bisimulation_test(fcrn, icrn, fs, interpretation = i)
        self.assertTrue(v)

    def todo_test_cardelliNM_modularity(self):
        # echo "->A; B->" | nuskell --ts cardelli2011_NM.ts --verify crn-bisimulation
        pass


if __name__ == '__main__':
    unittest.main()
