#!/usr/bin/env python
#
#  tests/test_crn_bisimulation.py
#  Original source from the Nuskell compiler project
#
import logging
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

import unittest

from collections import Counter
from crnverifier.utils import parse_crn
from crnverifier.crn_bisimulation import (crn_bisimulation_test, 
                                          modular_crn_bisimulation_test,
                                          crn_bisimulations)
from crnverifier.crn_bisimulation import (is_modular, 
                                          same_reaction, 
                                          enum,
                                          updateT, checkT)
from crnverifier.deprecated import moduleCond

SKIP_SLOW = True

class JustCuriousTests(unittest.TestCase):
    # Maybe those tests are useless, but I wonder ...
    def test_me_quickly_01(self):
        fcrn = "A + B -> C"
        icrn = "x + y -> c + d"
        fcrn, _ = parse_crn(fcrn)
        icrn, _ = parse_crn(icrn)

        print()
        #for e, b in enumerate(crn_bisimulations(fcrn, icrn), 1):
        #    print(e, b)

        assert len(list(crn_bisimulations(fcrn, icrn))) == 4

        fcrn = "A + B -> C"
        icrn = "x + y + z -> c + d"
        fcrn, _ = parse_crn(fcrn)
        icrn, _ = parse_crn(icrn)
        for e, b in enumerate(crn_bisimulations(fcrn, icrn), 1):
            #print(e, b)
            assert b is None

        fcrn = " -> A"
        icrn = " y -> a"
        fcrn, _ = parse_crn(fcrn)
        icrn, _ = parse_crn(icrn)
        for e, b in enumerate(crn_bisimulations(fcrn, icrn), 1):
            #print(e, b)
            assert b is None

        fcrn = " -> A"
        icrn = " -> y; y -> a"
        fcrn, _ = parse_crn(fcrn)
        icrn, _ = parse_crn(icrn)
        for e, b in enumerate(crn_bisimulations(fcrn, icrn), 1):
            #print(e, b)
            self.assertDictEqual(b, {'y': ['A'], 'a': ['A']})

        # TODO: this returns the same interpretation 3 times!
        fcrn = " -> A"
        icrn = " -> y; y <=> z; z -> a"
        fcrn, _ = parse_crn(fcrn)
        icrn, _ = parse_crn(icrn)
        for e, b in enumerate(crn_bisimulations(fcrn, icrn), 1):
            print(e, b)

        # TODO: this returns both interpretations 3 times!
        fcrn = "A -> "
        icrn = "a -> y; y <=> z; z -> "
        fcrn, _ = parse_crn(fcrn)
        icrn, _ = parse_crn(icrn)
        for e, b in enumerate(crn_bisimulations(fcrn, icrn), 1):
            print(e, b)


class HelperTests(unittest.TestCase):
    def test_enum(self):

        assert list(enum(0, Counter())) == [[]]
        assert list(enum(0, Counter(), [])) == [[]]
        assert list(enum(1, Counter(['C']), [1])) == [[Counter({'C': 1})]]
        assert list(enum(1, Counter(['a', 'b']))) == [[Counter({'a': 1, 'b': 1})]]

        out = [[Counter(), Counter({'B': 1, 'D': 1})], 
               [Counter({'D': 1}), Counter({'B': 1})],
               [Counter({'B': 1}), Counter({'D': 1})], 
               [Counter({'B': 1, 'D': 1}), Counter()]]
        assert list(enum(2, Counter(['B', 'D']), [1, 1])) == out

        #print(list(enum(3, Counter(['A', 'B', 'C']), [1, 1, 1])))

        with self.assertRaises(IndexError):
            # one weight, two parts.
            list(enum(2, Counter(['a']), [1]))


    def test_same_reaction(self):
        # literally same reaction.
        frxn = "A + B -> C"
        irxn = "A + B -> C"
        fcrn, fs = parse_crn(frxn)
        icrn, _ = parse_crn(irxn)
        frxn = [Counter(part) for part in fcrn[0]]
        irxn = [Counter(part) for part in icrn[0]]
        assert same_reaction(irxn, frxn, fs) is True

        # trying to break the old code ... 
        frxn = "A + B -> C + B"
        irxn = "A + b -> B"
        fcrn, fs = parse_crn(frxn)
        icrn, _ = parse_crn(irxn)
        frxn = [Counter(part) for part in fcrn[0]]
        irxn = [Counter(part) for part in icrn[0]]
        assert same_reaction(irxn, frxn, fs) is False

    def test_same_reaction_wrong(self):
        # trying to break the old code ... 
        frxn = "A -> C + D"
        irxn = "A + y -> C + y"
        fcrn, fs = parse_crn(frxn)
        icrn, _ = parse_crn(irxn)
        frxn = [Counter(part) for part in fcrn[0]]
        irxn = [Counter(part) for part in icrn[0]]
        assert same_reaction(irxn, frxn, fs) is False

        # trying to break the old code ... 
        frxn = "A -> C"
        irxn = "A + y -> C + y"
        fcrn, fs = parse_crn(frxn)
        icrn, _ = parse_crn(irxn)
        frxn = [Counter(part) for part in fcrn[0]]
        irxn = [Counter(part) for part in icrn[0]]
        assert same_reaction(irxn, frxn, fs) is True

    def test_same_reaction_products(self):
        # product interpretation 1
        frxn = "A + B -> C + D"
        irxn = "A + B -> c"
        fcrn, fs = parse_crn(frxn)
        icrn, _ = parse_crn(irxn)
        frxn = [Counter(part) for part in fcrn[0]]
        irxn = [Counter(part) for part in icrn[0]]
        assert same_reaction(irxn, frxn, fs) is True

        # product interpretation 2
        frxn = "A + B -> C"
        irxn = "A + B -> c + d"
        fcrn, fs = parse_crn(frxn)
        icrn, _ = parse_crn(irxn)
        frxn = [Counter(part) for part in fcrn[0]]
        irxn = [Counter(part) for part in icrn[0]]
        assert same_reaction(irxn, frxn, fs) is True

        # product interpretation 3
        frxn = "A + B -> C"
        irxn = "A + B -> "
        fcrn, fs = parse_crn(frxn)
        icrn, _ = parse_crn(irxn)
        frxn = [Counter(part) for part in fcrn[0]]
        irxn = [Counter(part) for part in icrn[0]]
        assert same_reaction(irxn, frxn, fs) is False

        # product interpretation 4
        frxn = "A + B -> C"
        irxn = "A + B -> C + B"
        fcrn, fs = parse_crn(frxn)
        icrn, _ = parse_crn(irxn)
        frxn = [Counter(part) for part in fcrn[0]]
        irxn = [Counter(part) for part in icrn[0]]
        assert same_reaction(irxn, frxn, fs) is False

    def test_same_reaction_reactants(self):
        frxn = "A + B -> C"
        irxn = "a -> C"
        fcrn, fs = parse_crn(frxn)
        icrn, _ = parse_crn(irxn)
        frxn = [Counter(part) for part in fcrn[0]]
        irxn = [Counter(part) for part in icrn[0]]
        assert same_reaction(irxn, frxn, fs) is True

        # Potential null species ...
        frxn = "A -> C + D"
        irxn = "a + b -> C + D"
        fcrn, fs = parse_crn(frxn)
        icrn, _ = parse_crn(irxn)
        frxn = [Counter(part) for part in fcrn[0]]
        irxn = [Counter(part) for part in icrn[0]]
        assert same_reaction(irxn, frxn, fs) is True

        # Potential null species ...
        frxn = "A -> C + D"
        irxn = "A + b -> C + D"
        fcrn, fs = parse_crn(frxn)
        icrn, _ = parse_crn(irxn)
        frxn = [Counter(part) for part in fcrn[0]]
        irxn = [Counter(part) for part in icrn[0]]
        assert same_reaction(irxn, frxn, fs) is True

        frxn = "A + B -> C"
        irxn = "A + B + A -> C"
        fcrn, fs = parse_crn(frxn)
        icrn, _ = parse_crn(irxn)
        frxn = [Counter(part) for part in fcrn[0]]
        irxn = [Counter(part) for part in icrn[0]]
        assert same_reaction(irxn, frxn, fs) is False

        frxn = "A + B -> C"
        irxn = "A -> C"
        fcrn, fs = parse_crn(frxn)
        icrn, _ = parse_crn(irxn)
        frxn = [Counter(part) for part in fcrn[0]]
        irxn = [Counter(part) for part in icrn[0]]
        assert same_reaction(irxn, frxn, fs) is False

    def test_update_table_01(self):
        fcrn = "A + B -> C"
        icrn = "x + y -> c + d"
        fcrn, fs = parse_crn(fcrn)
        icrn, _ = parse_crn(icrn)
        fcrn = [[Counter(part) for part in rxn] for rxn in fcrn]
        icrn = [[Counter(part) for part in rxn] for rxn in icrn]

        assert updateT(fcrn, icrn, fs) == [[True, True]]

        fcrn = "A + B -> C"
        icrn = "A + B -> C + d"
        fcrn, fs = parse_crn(fcrn)
        icrn, _ = parse_crn(icrn)
        fcrn = [[Counter(part) for part in rxn] for rxn in fcrn]
        icrn = [[Counter(part) for part in rxn] for rxn in icrn]
        assert updateT(fcrn, icrn, fs) == [[True, False]] 

        # TODO: this returns the same interpretation 3 times!
        fcrn = " -> A"
        icrn = " -> y; y <=> z; z -> a"
        fcrn, _ = parse_crn(fcrn)
        icrn, _ = parse_crn(icrn)
        fcrn = [[Counter(part) for part in rxn] for rxn in fcrn]
        icrn = [[Counter(part) for part in rxn] for rxn in icrn]
        assert updateT(fcrn, icrn, fs) == [[True, True], 
                                          [True, True],
                                          [True, True], 
                                          [True, True]]

        # TODO: this returns the same interpretation 3 times!
        fcrn = " -> A"
        icrn = " -> A; A <=> A; A -> A"
        fcrn, _ = parse_crn(fcrn)
        icrn, _ = parse_crn(icrn)
        fcrn = [[Counter(part) for part in rxn] for rxn in fcrn]
        icrn = [[Counter(part) for part in rxn] for rxn in icrn]
        assert updateT(fcrn, icrn, fs) == [[True, False],
                                          [False, True], 
                                          [False, True],
                                          [False, True]]

        # TODO: this returns the same interpretation 3 times!
        fcrn = " -> A"
        icrn = " -> ;  <=> ;  -> A"
        fcrn, _ = parse_crn(fcrn)
        icrn, _ = parse_crn(icrn)
        fcrn = [[Counter(part) for part in rxn] for rxn in fcrn]
        icrn = [[Counter(part) for part in rxn] for rxn in icrn]
        assert updateT(fcrn, icrn, fs) == [[False, True],
                                          [False, True],
                                          [False, True],
                                          [True, False]]

    def test_update_table_large(self):
        fcrn = [[Counter({'A': 1, 'b': 1}), Counter({'c': 1})], 
                [Counter({'b': 1}), Counter({'c': 1})],
                [Counter({'c': 1}), Counter({'b': 1})],
                [Counter({'b': 1}), Counter({'b': 2})]]
        icrn = [[Counter({'A': 1}), Counter({'i7': 1})],
                [Counter({'i7': 1}), Counter({'A': 1})],
                [Counter({'i7': 1, 'b': 1}), Counter({'i19': 1})],
                [Counter({'b': 1}), Counter({'i96': 1})],
                [Counter({'b': 1}), Counter({'i148': 1})],
                [Counter({'i7': 1, 'b': 1}), Counter({'i19': 1})],
                [Counter({'b': 1}), Counter({'i96': 1})],
                [Counter({'b': 1}), Counter({'i148': 1})],
                [Counter({'i7': 1, 'b': 1}), Counter({'i19': 1})],
                [Counter({'b': 1}), Counter({'i96': 1})],
                [Counter({'b': 1}), Counter({'i148': 1})],
                [Counter({'c': 1}), Counter({'i340': 1})], 
                [Counter({'c': 1}), Counter({'i340': 1})], 
                [Counter({'i19': 1}), Counter({'c': 1})],
                [Counter({'i96': 1}), Counter({'c': 1})],
                [Counter({'i148': 1}), Counter({'b': 2})],
                [Counter({'i340': 1}), Counter({'b': 1})]]
        fs = {'c', 'b', 'A'}
        result = [[False, False, False, False, True],
                  [False, False, False, False, True],
                  [True,  True,  False, True,  True], 
                  [False, True,  False, True,  True],
                  [False, True,  False, True,  True],
                  [True,  True,  False, True,  True], 
                  [False, True,  False, True,  True], 
                  [False, True,  False, True,  True],
                  [True,  True,  False, True,  True],
                  [False, True,  False, True,  True], 
                  [False, True,  False, True,  True],
                  [False, False, True,  False, True],
                  [False, False, True,  False, True], 
                  [True,  True,  False, False, True],
                  [True,  True,  False, False, True], 
                  [False, False, False, True,  True],
                  [False, False, True,  False, True]]
        assert updateT(fcrn, icrn, fs) == result

    def test_check_table(self):
        table = [[True, True]]
        assert checkT(table) is True
        table = [[False, False]]
        assert checkT(table) is False
        table = [[True, True],
                 [False, False]]
        assert checkT(table) is False
        table = [[False, True],
                 [True, False]]
        assert checkT(table) is True
        table = [[False, False],
                 [True, True]]
        assert checkT(table) is False
        table = [[True, False],
                 [True, False]]
        assert checkT(table) is True
        table = [[True, True, False],
                 [False, True, False]]
        assert checkT(table) is True

    def test_modularity_example_01(self):
        module = """ a <=> i1
                     b + i1 -> i2 + w3
                     i2 -> c + w4
                 """
        module, _ = parse_crn(module)
        fsc = {'A', 'B', 'C'}
        isc = {'a', 'b', 'c'}

        bisim = {'a': ['A'], 'b': ['B'], 'c': ['C'], 'i1': ['A'], 'i2': ['C'], 'w3': [], 'w4': []}
        assert moduleCond(module, fsc, isc, bisim) is True
        assert is_modular(bisim, module, isc, fsc) is True
        bisim = {'a': ['B'], 'b': ['A'], 'c': ['C'], 'i1': ['B'], 'i2': ['C'], 'w3': [], 'w4': []}
        assert moduleCond(module, fsc, isc, bisim) is True
        assert is_modular(bisim, module, isc, fsc) is True
        bisim = {'a': ['A'], 'b': ['B'], 'c': ['C'], 'i1': ['A'], 'i2': ['A','B'], 'w3': [], 'w4': []}
        assert moduleCond(module, fsc, isc, bisim) is False
        assert is_modular(bisim, module, isc, fsc) is False

    def test_modularity_example_02(self):
        module = """ b <=> e1
                     e1 -> e2 + e3 + e4
                     e4 -> e1 + e5
                 """
        module, _ = parse_crn(module)
        fsc = {'B'}
        isc = {'b'}

        minter = {'b': ['B'], 'e1': ['B'], 'e2': [], 'e3': [], 'e4': [], 'e5': []} 
        assert moduleCond(module, fsc, isc, minter) is True
        assert is_modular(minter, module, isc, fsc) is True

        minter = {'b': ['B'], 'e1': ['B'], 'e2': [], 'e3': [], 'e4': ['B'], 'e5': []} 
        assert moduleCond(module, fsc, isc, minter) is True
        assert is_modular(minter, module, isc, fsc) is True

        minter = {'b': ['B'], 'e1': ['B', 'B'], 'e2': [], 'e3': [], 'e4': [], 'e5': []} 
        assert moduleCond(module, fsc, isc, minter) is False
        assert is_modular(minter, module, isc, fsc) is False

        minter = {'b': ['B'], 'e1': ['B'], 'e2': ['B'], 'e3': [], 'e4': [], 'e5': []} 
        assert moduleCond(module, fsc, isc, minter) is False
        assert is_modular(minter, module, isc, fsc) is False
        isc = {'b', 'e2'}
        assert moduleCond(module, fsc, isc, minter) is True
        assert is_modular(minter, module, isc, fsc) is True
        isc = {'b', 'e5'}
        minter = {'b': ['B'], 'e1': [], 'e2': [], 'e3': [], 'e4': [], 'e5': ['B']} 
        assert moduleCond(module, fsc, isc, minter) is True
        assert is_modular(minter, module, isc, fsc) is True

        minter = {'b': ['B'], 'e1': ['B'], 'e2': ['B'], 'e3': [], 'e4': [], 'e5': ['B']} 
        assert moduleCond(module, fsc, isc, minter) is False
        assert is_modular(minter, module, isc, fsc) is False

    def test_modularity_example_03(self):
        module = """
                    a <=> e6
                    a <=> e1
                    a + e6 <=> e1
                """
        module, _ = parse_crn(module)
        fsc, isc = {'A'}, {'a'}
        inter = {'a': ['A'], 'e1': ['A', 'A'], 'e6': ['A']}
        assert moduleCond(module, fsc, isc, inter) is True
        assert is_modular(inter, module, isc, fsc) is True

    def test_modularity_example_04(self):
        fcrn = "A + B -> C + D"
        icrn = "a + b -> x; x -> y + z; y -> c; z -> d"
        fcrn, _ = parse_crn(fcrn)
        icrn, _ = parse_crn(icrn)

        print() 
        fsc, isc = {'A','B', 'C', 'D'}, {'a', 'b', 'c', 'd'}
        for e, bisim in enumerate(crn_bisimulations(fcrn, icrn)):
            print(e, bisim)
            assert moduleCond(icrn, fsc, isc, bisim) == is_modular(bisim, icrn, isc, fsc)

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

        #assert all(i in inters for i in [i01, i02, i03, i04, i05, i06, i07, i08, i09, i10])


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
