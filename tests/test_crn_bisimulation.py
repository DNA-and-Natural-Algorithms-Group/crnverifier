#!/usr/bin/env python
#
#  tests/test_crn_bisimulation.py
#  Original source from the Nuskell compiler project
#
import logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)

import unittest

from collections import Counter
from crnverifier.utils import parse_crn
from crnverifier.crn_bisimulation import (SpeciesAssignmentError,
                                          # Main interface
                                          crn_bisimulations,
                                          crn_bisimulation_test, 
                                          modular_crn_bisimulation_test,
                                          # HelperTests
                                          minimal_implementation_states,
                                          subsetsL, 
                                          enumL, 
                                          same_reaction, 
                                          updateT, 
                                          checkT,
                                          # ConditionTests
                                          passes_atomic_condition,
                                          passes_delimiting_condition,
                                          passes_permissive_condition,
                                          # Individual test classes
                                          search_column, 
                                          search_row, 
                                          passes_modularity_condition, 
                                          # Just used
                                          inter_counter,
                                          inter_list,
                                          subst) 

SKIP_SLOW = True
SKIP_DEBUG = False

def rl(rxn):
    return [list(part.elements()) for part in rxn]

def rc(rxn):
    return [Counter(part) for part in rxn]


@unittest.skipIf(SKIP_DEBUG, "skipping tests for debugging")
class JustCuriousTests(unittest.TestCase):
    # Some small examples that are easy to verify.
    def test_me_quickly_01(self):
        fcrn = "A + B -> C"
        icrn = "x + y -> c + d"
        fcrn, _ = parse_crn(fcrn)
        icrn, _ = parse_crn(icrn)
        bisims = list(crn_bisimulations(fcrn, icrn))
        if len(bisims) != 4:
            print('FAILURE:')
            for e, b in enumerate(bisims, 1):
                print(e, b)
        assert len(bisims) == 4
        assert bisims == list(crn_bisimulations(fcrn, icrn, permissive = 'graphsearch'))
        assert bisims == list(crn_bisimulations(fcrn, icrn, permissive = 'loopsearch'))
        assert bisims == list(crn_bisimulations(fcrn, icrn, permissive = 'bruteforce'))

    def test_me_quickly_02(self):
        fcrn = " -> A"
        icrn = " -> y; y -> a"
        fcrn, _ = parse_crn(fcrn)
        icrn, _ = parse_crn(icrn)
        bisims = list(crn_bisimulations(fcrn, icrn))
        b1 = {'a': ['A'], 'y': ['A']}
        b2 = {'a': ['A'], 'y': []}
        assert len(bisims) == 2
        assert b1 in bisims
        assert b2 in bisims
        assert bisims == list(crn_bisimulations(fcrn, icrn, permissive = 'graphsearch'))
        assert bisims == list(crn_bisimulations(fcrn, icrn, permissive = 'loopsearch'))
        assert bisims == list(crn_bisimulations(fcrn, icrn, permissive = 'bruteforce'))

    def test_me_quickly_03(self):
        fcrn = " -> A"
        icrn = " -> y; y <=> z; z -> a"
        fcrn, _ = parse_crn(fcrn)
        icrn, _ = parse_crn(icrn)
        bisims = list(crn_bisimulations(fcrn, icrn))
        b1 = {'a': ['A'], 'y': ['A'], 'z': ['A']}
        b2 = {'a': ['A'], 'y': [], 'z': []}
        if len(bisims) != 2:
            print('FAILURE:')
            for e, b in enumerate(bisims):
                print(e, b)
        assert len(bisims) == 2
        assert b1 in bisims
        assert b2 in bisims
        assert bisims == list(crn_bisimulations(fcrn, icrn, permissive = 'graphsearch'))
        assert bisims == list(crn_bisimulations(fcrn, icrn, permissive = 'loopsearch'))
        assert bisims == list(crn_bisimulations(fcrn, icrn, permissive = 'bruteforce'))

    def test_me_quickly_04(self):
        fcrn = "A -> "
        icrn = "a -> y; y <=> z; z -> "
        fcrn, _ = parse_crn(fcrn)
        icrn, _ = parse_crn(icrn)
        bisims = list(crn_bisimulations(fcrn, icrn))
        if len(bisims) != 2:
            print('FAILURE:')
            for e, b in enumerate(bisims):
                print(e, b)
        b1 = {'a': ['A'], 'y': [], 'z': []}
        b2 = {'a': ['A'], 'y': ['A'], 'z': ['A']}
        assert len(bisims) == 2
        assert b1 in bisims
        assert b2 in bisims
        assert bisims == list(crn_bisimulations(fcrn, icrn, permissive = 'graphsearch'))
        assert bisims == list(crn_bisimulations(fcrn, icrn, permissive = 'loopsearch'))
        assert bisims == list(crn_bisimulations(fcrn, icrn, permissive = 'bruteforce'))

    def test_me_quickly_false(self):
        fcrn = "A + B -> C"
        icrn = "x + y + z -> c + d"
        fcrn, _ = parse_crn(fcrn)
        icrn, _ = parse_crn(icrn)
        bisims = list(crn_bisimulations(fcrn, icrn))
        assert len(bisims) == 0
        assert bisims == list(crn_bisimulations(fcrn, icrn, permissive = 'graphsearch'))
        assert bisims == list(crn_bisimulations(fcrn, icrn, permissive = 'loopsearch'))
        assert bisims == list(crn_bisimulations(fcrn, icrn, permissive = 'bruteforce'))

        fcrn = " -> A"
        icrn = " y -> a"
        fcrn, _ = parse_crn(fcrn)
        icrn, _ = parse_crn(icrn)
        bisims = list(crn_bisimulations(fcrn, icrn))
        assert len(bisims) == 0
        assert bisims == list(crn_bisimulations(fcrn, icrn, permissive = 'graphsearch'))
        assert bisims == list(crn_bisimulations(fcrn, icrn, permissive = 'loopsearch'))
        assert bisims == list(crn_bisimulations(fcrn, icrn, permissive = 'bruteforce'))

@unittest.skipIf(SKIP_DEBUG, "skipping tests for debugging")
class HelperTests(unittest.TestCase):
    """ Helper functions for CRN bisimulation:
      - minimal_implementation_states
      - subsetsL
      - enumL
      - same_reaction
      - updateT
      - checkT
    """
    def test_minimal_implementation_states_exact(self):
        state = []
        inter = {'a': ['A'], 
                 'b': ['B'],
                 'x': ['A', 'B'],
                 'y': ['A', 'A']}
        assert list(minimal_implementation_states(state, inter)) == [[]]

        state = ['A']
        inter = {'a': ['A'], 
                 'b': ['B'],
                 'x': ['A', 'B'],
                 'y': ['A', 'A']}
        minis = [['a'], ['x'], ['y']]
        assert list(minimal_implementation_states(state, inter)) == minis

    def test_minimal_implementation_states_supersets(self):
        # NOTE: these used to return supersets, but now
        # that should be fixed.
        state = ['A', 'A']
        inter = {'a': ['A'], 
                 'b': ['B'],
                 'x': ['A', 'B'],
                 'y': ['A', 'A']}
        minis = [['a', 'a'], ['a', 'x'], ['x', 'x'], ['y']]
        supfs = list(map(sorted, minimal_implementation_states(state, inter))) 
        assert all(sorted(m) in supfs for m in minis)
        assert len(minis) == len(supfs)

        state = []
        inter = {'a': [], 
                 'b': []}
        minis = [[]]
        supfs = list(map(sorted, minimal_implementation_states(state, inter))) 
        assert all(sorted(m) in supfs for m in minis)
        assert len(minis) == len(supfs)


        state = ['A', 'B']
        inter = {'a': ['A'], 
                 'b': ['B'],
                 'c': ['C'],
                 'x': ['A', 'B']}
        minis = [['a', 'b'], ['x']]
        supfs = list(map(sorted, minimal_implementation_states(state, inter))) 
        assert all(sorted(m) in supfs for m in minis)
        assert len(minis) == len(supfs)

        state = ['A', 'B', 'A']
        inter = {'a': ['A'], 
                 'b': ['B'],
                 'c': ['C'],
                 'x': ['A', 'B']}
        minis = [['a', 'a', 'b'], ['x', 'a'], ['x', 'x']]
        supfs = list(map(sorted, minimal_implementation_states(state, inter))) 
        assert all(sorted(m) in supfs for m in minis)
        assert len(minis) == len(supfs)

    def test_subsets(self):
        #['A']      => [[], ['A']]
        #['A', 'A'] => [[], ['A'], ['A', 'A']
        #['A', 'B'] => [[], ['A'], ['B'], ['A', 'B']
        assert sorted(subsetsL([])) == [()]
        assert sorted(subsetsL(['A'])) == [(), ('A',)]
        assert sorted(subsetsL(['A', 'A'])) == [(), ('A',), ('A',), ('A', 'A')]
        assert sorted(set(subsetsL(['A', 'A']))) == [(), ('A',), ('A', 'A')]
        assert sorted(subsetsL(['A', 'B'])) == sorted([(), ('A',), ('B',), ('A', 'B')])
        assert sorted(set(subsetsL(['A', 'A', 'B']))) == sorted(
                [(), ('A',), ('B',), 
                 ('A', 'A'), ('A', 'B'), 
                 ('A', 'A', 'B')])

    def test_enum_noweights(self):
        # For example: 
        #   - n = 3 for three unassinged implmentation species (x, y, z).
        assert list(enumL(0, [])) == [[]]
        assert list(enumL(1, [])) == [[()]]
        assert list(enumL(2, [])) == [[(), ()]]
        assert list(enumL(3, [])) == [[(), (), ()]]
        assert list(enumL(4, [])) == [[(), (), (), ()]]

        assert list(enumL(0, ['A'])) == [[]]
        assert list(enumL(1, ['A'])) == [[('A',)]]
        assert list(enumL(1, ['A', 'B', 'C'])) == [[('A', 'B', 'C')]]
        assert sorted(enumL(2, ['A'])) == sorted([[('A',), ()], [(), ('A',)]])
        assert sorted(enumL(3, ['A', 'B'])) == sorted([
                 [('A',), ('B',), ()],
                 [('A',), (), ('B',)],
                 [('B',), ('A',), ()], 
                 [('B',), (), ('A',)],
                 [(), ('A',), ('B',)], 
                 [(), ('B',), ('A',)], 
                 [('A', 'B'), (), ()],
                 [(), ('A', 'B'), ()],
                 [(), (), ('A', 'B')]])

        assert sorted(enumL(2, ['A', 'A', 'B'])) == sorted([ 
                 [(), ('A', 'A', 'B')], 
                 [('A', 'A', 'B'), ()], 
                 [('A', 'B'), ('A',)],
                 [('A', 'A'), ('B',)], 
                 [('B',), ('A', 'A')], 
                 [('A',), ('A', 'B')]])

        assert sorted(enumL(2, ['A', 'B', 'C'])) == sorted([
                 [('A', 'B', 'C'), ()], 
                 [(), ('A', 'B', 'C')], 
                 [('A',), ('B', 'C')], 
                 [('B',), ('A', 'C')], 
                 [('C',), ('A', 'B')], 
                 [('A', 'B'), ('C',)], 
                 [('A', 'C'), ('B',)],
                 [('B', 'C'), ('A',)]])

    def test_enum_weights(self):
        assert list(enumL(0, [], weights = [])) == [[]]
        assert list(enumL(1, [], weights = [1])) == [[()]]
        assert list(enumL(2, [], weights = [1, 1])) == [[(), ()]]
        assert list(enumL(3, [], weights = [1, 2, 3])) == [[(), (), ()]]

        assert list(enumL(1, ['A'], weights = [1])) == [[('A',)]]
        with self.assertRaises(SpeciesAssignmentError):
            assert list(enumL(1, ['A'], weights = [2])) == [[()]]
        assert list(enumL(1, ['A', 'A'], weights = [2])) == [[('A',)]]
        assert sorted(list(enumL(2, ['A', 'A'], weights = [2, 1]))) == sorted(
                [[('A',), ()], [(), ('A', 'A')]])
        assert sorted(list(enumL(2, ['A', 'A'], weights = [1, 2]))) == sorted(
                [[(), ('A',)], [('A', 'A'), ()]])
 
        assert sorted(list(enumL(2, ['A', 'B'], weights = [1, 2]))) == sorted(
                [[('A', 'B'), ()]])

        assert sorted(list(enumL(2, ['A', 'A', 'B'], weights = [1, 2]))) == sorted(
                [[('A', 'A', 'B'), ()],
                 [('B',), ('A',)]])

        assert sorted(list(enumL(3, list('AAAAB'), weights = [2, 1, 2]))) == sorted([
                    [(), ('A', 'A', 'A', 'A', 'B'), ()], 
                    [('A',), ('A', 'A', 'B'), ()],
                    [(), ('A', 'A', 'B'), ('A',)], 
                    [('A', 'A'), ('B',), ()],
                    [(), ('B',), ('A', 'A')],
                    [('A',), ('B',), ('A',)]])

    def test_same_reaction(self):
        frxn = "A + B -> C"
        irxn = "A + B -> C"
        fcrn, fs = parse_crn(frxn)
        icrn, _ = parse_crn(irxn)
        assert same_reaction(icrn[0], fcrn[0], fs, counter = False)

        frxn = "A + B -> C + B"
        irxn = "A + b -> B"
        fcrn, fs = parse_crn(frxn)
        icrn, _ = parse_crn(irxn)
        assert not same_reaction(icrn[0], fcrn[0], fs, counter = False)

    def test_same_reaction_new(self):
        # trying to break the old code ... 
        frxn = "A -> C + D"
        irxn = "A + y -> C + y"
        fcrn, fs = parse_crn(frxn)
        icrn, _ = parse_crn(irxn)
        assert not same_reaction(icrn[0], fcrn[0], fs, counter = False)

        frxn = "A -> C"
        irxn = "A + y -> C + y"
        fcrn, fs = parse_crn(frxn)
        icrn, _ = parse_crn(irxn)
        assert same_reaction(icrn[0], fcrn[0], fs, counter = False)

    def test_same_reaction_products(self):
        frxn = "A + B -> C + D"
        irxn = "A + B -> c"
        fcrn, fs = parse_crn(frxn)
        icrn, _ = parse_crn(irxn)
        assert same_reaction(icrn[0], fcrn[0], fs, counter = False)

        frxn = "A + B -> C"
        irxn = "A + B -> c + d"
        fcrn, fs = parse_crn(frxn)
        icrn, _ = parse_crn(irxn)
        assert same_reaction(icrn[0], fcrn[0], fs, counter = False)

        frxn = "A + B -> C"
        irxn = "A + B -> "
        fcrn, fs = parse_crn(frxn)
        icrn, _ = parse_crn(irxn)
        assert not same_reaction(icrn[0], fcrn[0], fs, counter = False)

        frxn = "A + B -> C"
        irxn = "A + B -> C + B"
        fcrn, fs = parse_crn(frxn)
        icrn, _ = parse_crn(irxn)
        assert not same_reaction(icrn[0], fcrn[0], fs, counter = False)

    def test_same_reaction_reactants(self):
        # NOTE: tests include potential null species ...
        frxn = "A + B -> C"
        irxn = "a -> C"
        fcrn, fs = parse_crn(frxn)
        icrn, _ = parse_crn(irxn)
        assert same_reaction(icrn[0], fcrn[0], fs, counter = False)

        frxn = "A -> C + D"
        irxn = "a + b -> C + D"
        fcrn, fs = parse_crn(frxn)
        icrn, _ = parse_crn(irxn)
        assert same_reaction(icrn[0], fcrn[0], fs, counter = False)

        frxn = "A -> C + D"
        irxn = "A + b -> C + D"
        fcrn, fs = parse_crn(frxn)
        icrn, _ = parse_crn(irxn)
        assert same_reaction(icrn[0], fcrn[0], fs, counter = False)

        frxn = "A + B -> C"
        irxn = "A + B + A -> C"
        fcrn, fs = parse_crn(frxn)
        icrn, _ = parse_crn(irxn)
        assert not same_reaction(icrn[0], fcrn[0], fs, counter = False)

        frxn = "A + B -> C"
        irxn = "A -> C"
        fcrn, fs = parse_crn(frxn)
        icrn, _ = parse_crn(irxn)
        assert not same_reaction(icrn[0], fcrn[0], fs, counter = False)

    def test_update_table(self):
        fcrn = "A + B -> C"
        icrn = "x + y -> c + d"
        fcrn, fs = parse_crn(fcrn)
        icrn, _ = parse_crn(icrn)
        table = [[True, True]]
        assert updateT(fcrn, icrn, fs, counter = False) == table

        fcrn = "A + B -> C"
        icrn = "A + B -> C + d"
        fcrn, fs = parse_crn(fcrn)
        icrn, _ = parse_crn(icrn)
        table = [[True, False]]
        assert updateT(fcrn, icrn, fs, counter = False) == table

        fcrn = " -> A"
        icrn = " -> y; y <=> z; z -> a"
        fcrn, fs = parse_crn(fcrn)
        icrn, _ = parse_crn(icrn)
        table = [[True, True], 
                 [True, True],
                 [True, True], 
                 [True, True]]
        assert updateT(fcrn, icrn, fs, counter = False) == table

        fcrn = " -> A"
        icrn = " -> A; A <=> A; A -> A"
        fcrn, fs = parse_crn(fcrn)
        icrn, _ = parse_crn(icrn)
        table = [[True, False],
                 [False, True], 
                 [False, True],
                 [False, True]]
        assert updateT(fcrn, icrn, fs, counter = False) == table

        fcrn = " -> A"
        icrn = " -> ;  <=> ;  -> A"
        fcrn, _ = parse_crn(fcrn)
        icrn, _ = parse_crn(icrn)
        table = [[False, True],
                 [False, True],
                 [False, True],
                 [True, False]]
        assert updateT(fcrn, icrn, fs, counter = False) == table

    def test_update_table_large(self):
        fcrn = """ A + b -> c
                   b -> c
                   c -> b
                   b -> 2b """
        icrn = """ A -> i7
                   i7 -> A
                   i7 + b -> i19
                   b -> i96
                   b -> i148
                   i7 + b -> i19
                   b -> i96
                   b -> i148
                   i7 + b -> i19
                   b -> i96
                   b -> i148
                   c -> i340
                   c -> i340
                   i19 -> c
                   i96 -> c
                   i148 -> b + b
                   i340 -> b """
        fcrn, fs = parse_crn(fcrn)
        icrn, _ = parse_crn(icrn)
        table = [[False, False, False, False, True],
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
        assert updateT(fcrn, icrn, fs, counter = False) == table

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

@unittest.skipIf(SKIP_DEBUG, "skipping tests for debugging")
class ConditionTests(unittest.TestCase):
    def test_atomic_01(self):
        fs = set(['A', 'B', 'C'])
        inter = {'a' : ['B'],
                 'B' : ['A'],
                 'c' : ['C']}
        inter = inter_counter(inter)
        assert passes_atomic_condition(inter, fs)

    def test_delimiting_01(self):
        fcrn = "A + B -> C"
        icrn = "a + b -> c"
        fcrn, fs = parse_crn(fcrn)
        icrn, _ = parse_crn(icrn)
        fcrn = [rc(rxn) for rxn in fcrn]
        icrn = [rc(rxn) for rxn in icrn]
        inter = {'a' : ['B'],
                 'b' : ['A'],
                 'c' : ['C']}
        inter = inter_counter(inter)
        assert passes_delimiting_condition(fcrn, icrn, fs, inter)

    def test_permissive_01(self):
        fcrn = "A + B -> C"
        icrn = "a + b -> c"
        fcrn, fs = parse_crn(fcrn)
        icrn, _ = parse_crn(icrn)
        fcrn = [rc(rxn) for rxn in fcrn]
        icrn = [rc(rxn) for rxn in icrn]
        inter = {'a' : ['B'],
                 'b' : ['A'],
                 'c' : ['C']}
        inter = inter_counter(inter)
        passes, info = passes_permissive_condition(fcrn, icrn, fs, inter, 'graphsearch')
        assert passes
        passes, info = passes_permissive_condition(fcrn, icrn, fs, inter, 'loopsearch')
        assert passes
        passes, info = passes_permissive_condition(fcrn, icrn, fs, inter, 'bruteforce')
        assert passes

    def test_permissive_02(self):
        fcrn = " -> A"
        icrn = "x -> a"
        fcrn, fs = parse_crn(fcrn)
        icrn, _ = parse_crn(icrn)
        fcrn = [rc(rxn) for rxn in fcrn]
        icrn = [rc(rxn) for rxn in icrn]
        inter = {'a' : ['A'],
                 'x' : []}
        inter = inter_counter(inter)
        passes, info = passes_permissive_condition(fcrn, icrn, fs, inter, 'graphsearch')
        assert not passes
        passes, info = passes_permissive_condition(fcrn, icrn, fs, inter, 'loopsearch')
        assert not passes
        passes, info = passes_permissive_condition(fcrn, icrn, fs, inter, 'bruteforce')
        assert not passes

    def test_permissive_03(self):
        fcrn = " -> A"
        icrn = " -> x; x -> a"
        fcrn, fs = parse_crn(fcrn)
        icrn, _ = parse_crn(icrn)
        fcrn = [rc(rxn) for rxn in fcrn]
        icrn = [rc(rxn) for rxn in icrn]
        inter = {'a' : ['A'],
                 'x' : ['A']}
        inter = inter_counter(inter)
        passes, info = passes_permissive_condition(fcrn, icrn, fs, inter, 'graphsearch')
        assert passes
        passes, info = passes_permissive_condition(fcrn, icrn, fs, inter, 'loopsearch')
        assert passes
        passes, info = passes_permissive_condition(fcrn, icrn, fs, inter, 'bruteforce')
        assert passes


        inter = {'a' : ['A'],
                 'x' : []}
        inter = inter_counter(inter)
        passes, info = passes_permissive_condition(fcrn, icrn, fs, inter, 'graphsearch')
        assert passes
        passes, info = passes_permissive_condition(fcrn, icrn, fs, inter, 'loopsearch')
        assert passes
        passes, info = passes_permissive_condition(fcrn, icrn, fs, inter, 'bruteforce')
        assert passes

    def test_permissive_03b(self):
        fcrn = " -> A"
        icrn = " -> x; -> y; x + y -> a"
        fcrn, fs = parse_crn(fcrn)
        icrn, _ = parse_crn(icrn)
        fcrn = [rc(rxn) for rxn in fcrn]
        icrn = [rc(rxn) for rxn in icrn]
        inter = {'a' : ['A'],
                 'x' : [],
                 'y' : []}
        inter = inter_counter(inter)
        passes, info = passes_permissive_condition(fcrn, icrn, fs, inter, 'graphsearch')
        assert passes
        passes, info = passes_permissive_condition(fcrn, icrn, fs, inter, 'loopsearch')
        assert passes
        passes, info = passes_permissive_condition(fcrn, icrn, fs, inter, 'bruteforce')
        assert passes

        inter = {'a' : ['A'],
                 'x' : ['A'],
                 'y' : []}
        inter = inter_counter(inter)
        passes, info = passes_permissive_condition(fcrn, icrn, fs, inter, 'graphsearch')
        assert passes
        passes, info = passes_permissive_condition(fcrn, icrn, fs, inter, 'loopsearch')
        assert passes
        passes, info = passes_permissive_condition(fcrn, icrn, fs, inter, 'bruteforce')
        assert passes

    def test_permissive_04(self):
        fcrn = "B -> A"
        icrn = "b -> b + x; b + 3x -> a"
        fcrn, fs = parse_crn(fcrn)
        icrn, _ = parse_crn(icrn)
        fcrn = [rc(rxn) for rxn in fcrn]
        icrn = [rc(rxn) for rxn in icrn]
        inter = {'a' : ['A'],
                 'b' : ['B'],
                 'x' : []}
        inter = inter_counter(inter)
        passes, info = passes_permissive_condition(fcrn, icrn, fs, inter, 'graphsearch')
        assert passes
        passes, info = passes_permissive_condition(fcrn, icrn, fs, inter, 'loopsearch')
        assert passes
        passes, info = passes_permissive_condition(fcrn, icrn, fs, inter, 'bruteforce')
        assert passes

    def test_permissive_05(self):
        fcrn = "A -> B"
        icrn = "x -> y; y -> x + z; x + 3z -> b"
        fcrn, fs = parse_crn(fcrn)
        icrn, _ = parse_crn(icrn)
        fcrn = [rc(rxn) for rxn in fcrn]
        icrn = [rc(rxn) for rxn in icrn]
        inter = {'x' : ['A'],
                 'b' : ['B'],
                 'y' : ['A'],
                 'z' : []}
        inter = inter_counter(inter)
        passes, info = passes_permissive_condition(fcrn, icrn, fs, inter, 'graphsearch')
        assert passes
        passes, info = passes_permissive_condition(fcrn, icrn, fs, inter, 'loopsearch')
        assert passes
        passes, info = passes_permissive_condition(fcrn, icrn, fs, inter, 'bruteforce')
        assert passes

    def test_permissive_06(self):
        fcrn = "A + B -> C + D"
        icrn = "a + b -> c + d; d + c -> e + f"
        fcrn, fs = parse_crn(fcrn)
        icrn, _ = parse_crn(icrn)
        fcrn = [rc(rxn) for rxn in fcrn]
        icrn = [rc(rxn) for rxn in icrn]
        inter={'a': ['B'], 'b': ['A'], 'c': [], 'd': ['A', 'B'], 'e': ['C'], 'f': ['D']}
        inter = inter_counter(inter)
        passes, info = passes_permissive_condition(fcrn, icrn, fs, inter, 'graphsearch')
        assert not passes
        passes, info = passes_permissive_condition(fcrn, icrn, fs, inter, 'loopsearch')
        assert not passes
        passes, info = passes_permissive_condition(fcrn, icrn, fs, inter, 'bruteforce')
        assert not passes

    def test_permissive_07(self):
        # JDW 2019
        fcrn = "A + B -> C"
        icrn = """ a1 <=> a2
                   a2 + b1 <=> ab
                   ab -> a1 + b1 + 2z
                   b1 + 3z -> b2
                   a1 + b2 + 2z -> c1
                   a2 + b2 -> c2
               """
        fcrn, fs = parse_crn(fcrn)
        icrn, _ = parse_crn(icrn)
        fcrn = [rc(rxn) for rxn in fcrn]
        icrn = [rc(rxn) for rxn in icrn]
        inter={'a1': ['A'], 
               'a2': ['A'], 
               'b1': ['B'],
               'b2': ['B'], 
               'ab': ['A', 'B'], 
               'c1': ['C'], 
               'c2': ['C'], 
               'z': []}
        inter = inter_counter(inter)
        passes, info = passes_permissive_condition(fcrn, icrn, fs, inter, 'graphsearch')
        assert passes
        passes, info = passes_permissive_condition(fcrn, icrn, fs, inter, 'loopsearch')
        assert passes
        passes, info = passes_permissive_condition(fcrn, icrn, fs, inter, 'bruteforce')
        assert passes

    def test_permissive_08(self):
        fcrn, fs = parse_crn('tests/crns/crn6.crn', is_file = True)
        icrn, _ = parse_crn('tests/crns/icrns/crn6_qingdong_thesis.crn', is_file = True)
        inter = {'i{A}': ['A'], 
                 'i{B}': ['B'], 
                 'i{C}': ['C'], 
                 'i{X}': ['X'],
                 'i{Y}': ['Y'], 
                 'i14': [],
                 'i15': [], 
                 'i73': ['B'], 
                 'i119': ['X', 'B', 'A'], 
                 'i120': [], 
                 'i194': [], 
                 'i394': ['X', 'X', 'Y'], 
                 'i575': ['X'], 
                 'i599': ['C'], 
                 'i631': [], 
                 'i778': ['Y'], 
                 'i842': ['Y', 'X', 'A'], 
                 'i886': [],
                 'i969': [],
                 'i1457': [], 
                 'i2232': ['A'], 
                 'i2300': ['A', 'C'], 
                 'i2340': [], 
                 'i2392': [], 
                 'i3032': []}
        fcrn = [rc(rxn) for rxn in fcrn]
        icrn = [rc(rxn) for rxn in icrn]
        inter = inter_counter(inter)
        passes, info = passes_permissive_condition(fcrn, icrn, fs, inter, 'graphsearch')
        assert not passes
        #passes, info = passes_permissive_condition(fcrn, icrn, fs, inter, 'loopsearch')
        #assert not passes
        passes, info = passes_permissive_condition(fcrn, icrn, fs, inter, 'bruteforce')
        assert not passes

@unittest.skipIf(SKIP_DEBUG, "skipping tests for debugging")
class TestColumnSearch(unittest.TestCase):
    def test_search_column_01(self):
        fcrn = "A -> B + C"
        icrn = """x1 -> x2
                  x2 -> x3 + x4
                  x3 <=> x5
                  x4 -> x7 + x8 
               """
        fcrn, fs = parse_crn(fcrn)
        icrn, _ = parse_crn(icrn)
        fcrn = [rc(rxn) for rxn in fcrn]
        icrn = [rc(rxn) for rxn in icrn]

        i1 = {'x2': ['A'], 'x3': ['B', 'C'], 'x4': []}
        i2 = {'x2': ['A'], 'x3': ['B'], 'x4': ['C']}
        i3 = {'x2': ['A'], 'x3': ['C'], 'x4': ['B']}
        i4 = {'x2': ['A'], 'x3': [], 'x4': ['B', 'C']}
        i5 = {'x4': ['A'], 'x7': ['B', 'C'], 'x8': []}
        i6 = {'x4': ['A'], 'x7': ['B'], 'x8': ['C']}
        i7 = {'x4': ['A'], 'x7': ['C'], 'x8': ['B']}
        i8 = {'x4': ['A'], 'x7': [], 'x8': ['B', 'C']}
        i9 = {'x1': ['A'], 'x2': ['B', 'C']}
        cols = list(search_column(fcrn, icrn, fs))
        if len(cols) != 9:
            print('FAILURE:')
            for e, b in enumerate(cols, 1):
                print(e, inter_list(b))
        assert len(cols) == 9
        assert inter_counter(i1) in cols
        assert inter_counter(i2) in cols
        assert inter_counter(i3) in cols
        assert inter_counter(i4) in cols
        assert inter_counter(i5) in cols
        assert inter_counter(i6) in cols
        assert inter_counter(i7) in cols
        assert inter_counter(i8) in cols
        assert inter_counter(i9) in cols

    def test_search_column_02(self):
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
        fcrn = [rc(rxn) for rxn in fcrn]
        icrn = [rc(rxn) for rxn in icrn]
        cols = list(search_column(fcrn, icrn, fs))

        # SB: I didn't actually check those, but you should if something goes wrong!
        i01 = {'x1': ['A', 'B'], 'x2': ['C', 'D'], 'x5': ['A', 'C'], 'x7': ['B', 'D']}
        i02 = {'x1': ['A', 'B'], 'x2': ['C', 'D'], 'x7': ['A', 'C'], 'x9': ['B', 'D'], 'x8': []}
        i03 = {'x1': ['A', 'B'], 'x2': ['C', 'D'], 'x7': ['A', 'C'], 'x9': ['B'], 'x8': ['D']}
        i04 = {'x2': ['A', 'B'], 'x6': ['C'], 'x8': ['D'], 'x5': ['A', 'C'], 'x7': ['B', 'D']}
        i05 = {'x2': ['A', 'B'], 'x6': ['C'], 'x8': ['D'], 'x7': ['A', 'C'], 'x9': ['B']}
        i06 = {'x2': ['A', 'B'], 'x6': ['C', 'D'], 'x8': [], 'x5': ['A', 'C'], 'x7': ['B', 'D']}
        i07 = {'x2': ['A', 'B'], 'x6': ['C', 'D'], 'x8': [], 'x7': ['A', 'C'], 'x9': ['B', 'D']}
        i08 = {'x5': ['A', 'B'], 'x7': ['C', 'D'], 'x1': ['A', 'C'], 'x2': ['B', 'D']}
        i09 = {'x5': ['A', 'B'], 'x7': ['C', 'D'], 'x2': ['A', 'C'], 'x6': ['B', 'D'], 'x8': []}
        i10 = {'x5': ['A', 'B'], 'x7': ['C', 'D'], 'x2': ['A', 'C'], 'x6': ['B'], 'x8': ['D']}
        i11 = {'x7': ['A', 'B'], 'x9': ['C'], 'x8': ['D'], 'x1': ['A', 'C'], 'x2': ['B', 'D']}
        i12 = {'x7': ['A', 'B'], 'x9': ['C'], 'x8': ['D'], 'x2': ['A', 'C'], 'x6': ['B']}
        i13 = {'x7': ['A', 'B'], 'x9': ['C', 'D'], 'x8': [], 'x1': ['A', 'C'], 'x2': ['B', 'D']}
        i14 = {'x7': ['A', 'B'], 'x9': ['C', 'D'], 'x8': [], 'x2': ['A', 'C'], 'x6': ['B', 'D']}
        i15 = {'x1': ['A', 'C'], 'x2': ['B', 'D'], 'x5': ['A', 'B'], 'x7': ['C', 'D']}
        i16 = {'x1': ['A', 'C'], 'x2': ['B', 'D'], 'x7': ['A', 'B'], 'x9': ['C'], 'x8': ['D']}
        i17 = {'x1': ['A', 'C'], 'x2': ['B', 'D'], 'x7': ['A', 'B'], 'x9': ['C', 'D'], 'x8': []}
        i18 = {'x2': ['A', 'C'], 'x6': ['B', 'D'], 'x8': [], 'x5': ['A', 'B'], 'x7': ['C', 'D']}
        i19 = {'x2': ['A', 'C'], 'x6': ['B', 'D'], 'x8': [], 'x7': ['A', 'B'], 'x9': ['C', 'D']}
        i20 = {'x2': ['A', 'C'], 'x6': ['B'], 'x8': ['D'], 'x5': ['A', 'B'], 'x7': ['C', 'D']}
        i21 = {'x2': ['A', 'C'], 'x6': ['B'], 'x8': ['D'], 'x7': ['A', 'B'], 'x9': ['C']}
        i22 = {'x5': ['A', 'C'], 'x7': ['B', 'D'], 'x1': ['A', 'B'], 'x2': ['C', 'D']}
        i23 = {'x5': ['A', 'C'], 'x7': ['B', 'D'], 'x2': ['A', 'B'], 'x6': ['C'], 'x8': ['D']}
        i24 = {'x5': ['A', 'C'], 'x7': ['B', 'D'], 'x2': ['A', 'B'], 'x6': ['C', 'D'], 'x8': []}
        i25 = {'x7': ['A', 'C'], 'x9': ['B', 'D'], 'x8': [], 'x1': ['A', 'B'], 'x2': ['C', 'D']}
        i26 = {'x7': ['A', 'C'], 'x9': ['B', 'D'], 'x8': [], 'x2': ['A', 'B'], 'x6': ['C', 'D']}
        i27 = {'x7': ['A', 'C'], 'x9': ['B'], 'x8': ['D'], 'x1': ['A', 'B'], 'x2': ['C', 'D']}
        i28 = {'x7': ['A', 'C'], 'x9': ['B'], 'x8': ['D'], 'x2': ['A', 'B'], 'x6': ['C']}

        if len(cols) != 28:
            print('FAILURE:')
            for e, b in enumerate(cols, 1):
                print(e, inter_list(b))
        assert len(cols) == 28 

@unittest.skipIf(SKIP_DEBUG, "skipping tests for debugging")
class TestRowSearch(unittest.TestCase):
    def test_search_row_01(self):
        fcrn = "A -> B + C"
        icrn = """x1 -> x2
                  x2 -> x3 + x4
                  x3 <=> x5
                  x4 -> x7 + x8 
               """
        fcrn, fs = parse_crn(fcrn)
        icrn, _ = parse_crn(icrn)
        fcrn = [rc(rxn) for rxn in fcrn]
        icrn = [rc(rxn) for rxn in icrn]
        i1 = {'x2': ['A'], 'x3': ['B', 'C'], 'x4': []}
        rows = list(search_row(fcrn, icrn, fs, inter_counter(i1)))
        assert len(rows) == 0

        i2 = {'x2': ['A'], 'x3': ['B'], 'x4': ['C']}
        i2r01 = {'x1': ['A'], 'x2': ['A'], 'x3': ['B'], 'x4': ['C'], 'x5': ['B'], 'x7': ['C'], 'x8': []}
        i2r02 = {'x1': ['A'], 'x2': ['A'], 'x3': ['B'], 'x4': ['C'], 'x5': ['B'], 'x7': [], 'x8': ['C']}
        rows = list(search_row(fcrn, icrn, fs, inter_counter(i2)))
        if len(rows) != 2:
            print('FAILURE:')
            for e, b in enumerate(rows, 1):
                print(e, {k: v for k, v in sorted(inter_list(b).items())})
        assert len(rows) == 2
        assert inter_counter(i2r01) in rows
        assert inter_counter(i2r02) in rows

        i3 = {'x2': ['A'], 'x3': ['C'], 'x4': ['B']}
        i3r01 = {'x1': ['A'], 'x2': ['A'], 'x3': ['C'], 'x4': ['B'], 'x5': ['C'], 'x7': [], 'x8': ['B']}
        i3r02 = {'x1': ['A'], 'x2': ['A'], 'x3': ['C'], 'x4': ['B'], 'x5': ['C'], 'x7': ['B'], 'x8': []}
        rows = list(search_row(fcrn, icrn, fs, inter_counter(i3)))
        if len(rows) != 2:
            print('FAILURE:')
            for e, b in enumerate(rows, 1):
                print(e, {k: v for k, v in sorted(inter_list(b).items())})
        assert len(rows) == 2
        assert inter_counter(i3r01) in rows
        assert inter_counter(i3r02) in rows

        i4 = {'x2': ['A'], 'x3': [], 'x4': ['B', 'C']}  
        i4r01 = {'x1': ['A'], 'x2': ['A'], 'x3': [], 'x4': ['B', 'C'], 'x5': [], 'x7': ['B'], 'x8': ['C']}
        i4r02 = {'x1': ['A'], 'x2': ['A'], 'x3': [], 'x4': ['B', 'C'], 'x5': [], 'x7': ['C'], 'x8': ['B']}
        rows = list(search_row(fcrn, icrn, fs, inter_counter(i4)))
        if len(rows) != 2:
            print('FAILURE:')
            for e, b in enumerate(rows, 1):
                print(e, {k: v for k, v in sorted(inter_list(b).items())})
        assert len(rows) == 2
        assert inter_counter(i4r01) in rows
        assert inter_counter(i4r02) in rows

        i5 = {'x4': ['A'], 'x7': ['B', 'C'], 'x8': []} 
        rows = list(search_row(fcrn, icrn, fs, inter_counter(i5)))
        assert len(rows) == 0

        i6 = {'x4': ['A'], 'x7': ['B'], 'x8': ['C']}
        i6r01 = {'x1': ['A'], 'x2': ['A'], 'x3': [], 'x4': ['A'], 'x5': [], 'x7': ['B'], 'x8': ['C']}
        rows = list(search_row(fcrn, icrn, fs, inter_counter(i6)))
        if len(rows) != 1:
            print('FAILURE:')
            for e, b in enumerate(rows, 1):
                print(e, {k: v for k, v in sorted(inter_list(b).items())})
        assert len(rows) == 1
        assert rows[0] == inter_counter(i6r01)
 
        i7 = {'x4': ['A'], 'x7': ['C'], 'x8': ['B']}
        i7r01 = {'x1': ['A'], 'x2': ['A'], 'x3': [], 'x4': ['A'], 'x5': [], 'x7': ['C'], 'x8': ['B']}
        rows = list(search_row(fcrn, icrn, fs, inter_counter(i7)))
        if len(rows) != 1:
            print('FAILURE:')
            for e, b in enumerate(rows, 1):
                print(e, {k: v for k, v in sorted(inter_list(b).items())})
        assert len(rows) == 1
        assert rows[0] == inter_counter(i7r01)

        i8 = {'x4': ['A'], 'x7': [], 'x8': ['B', 'C']} 
        rows = list(search_row(fcrn, icrn, fs, inter_counter(i8)))
        assert len(rows) == 0

        i9 = {'x1': ['A'], 'x2': ['B', 'C']} 
        i9r01 = {'x1': ['A'], 'x2': ['B', 'C'], 'x3': ['C'], 'x4': ['B'], 'x5': ['C'], 'x7': [], 'x8': ['B']}
        i9r02 = {'x1': ['A'], 'x2': ['B', 'C'], 'x3': ['C'], 'x4': ['B'], 'x5': ['C'], 'x7': ['B'], 'x8': []}
        i9r03 = {'x1': ['A'], 'x2': ['B', 'C'], 'x3': ['B'], 'x4': ['C'], 'x5': ['B'], 'x7': ['C'], 'x8': []}
        i9r04 = {'x1': ['A'], 'x2': ['B', 'C'], 'x3': ['B'], 'x4': ['C'], 'x5': ['B'], 'x7': [], 'x8': ['C']}
        i9r05 = {'x1': ['A'], 'x2': ['B', 'C'], 'x3': [], 'x4': ['B', 'C'], 'x5': [], 'x7': ['C'], 'x8': ['B']}
        i9r06 = {'x1': ['A'], 'x2': ['B', 'C'], 'x3': [], 'x4': ['B', 'C'], 'x5': [], 'x7': ['B'], 'x8': ['C']}
        rows = list(search_row(fcrn, icrn, fs, inter_counter(i9)))
        if len(rows) != 6:
            print('FAILURE:')
            for e, b in enumerate(rows, 1):
                print(e, {k: v for k, v in sorted(inter_list(b).items())})
        assert len(rows) == 6
        assert inter_counter(i9r01) in rows
        assert inter_counter(i9r02) in rows
        assert inter_counter(i9r03) in rows
        assert inter_counter(i9r04) in rows
        assert inter_counter(i9r05) in rows
        assert inter_counter(i9r06) in rows

@unittest.skipIf(SKIP_DEBUG, "skipping tests for debugging")
class TestSearchSpace(unittest.TestCase):
    def test_1f_1i(self):
        fcrn = "A + B -> C + D"
        icrn = "a + b -> c + d"
        fcrn, fs = parse_crn(fcrn)
        icrn, _ = parse_crn(icrn)

        i01 = {'a': ['A'], 'b': ['B'], 'c': ['C'], 'd': ['D']}
        i02 = {'a': ['B'], 'b': ['A'], 'c': ['C'], 'd': ['D']}
        i03 = {'a': ['A'], 'b': ['B'], 'c': ['D'], 'd': ['C']}
        i04 = {'a': ['B'], 'b': ['A'], 'c': ['D'], 'd': ['C']}
        bisims = list(crn_bisimulations(fcrn, icrn))
        assert len(bisims) == 4 
        assert i01 in bisims
        assert i02 in bisims
        assert i03 in bisims
        assert i04 in bisims

    def test_1f_2i(self):
        fcrn = "A + B -> C + D"
        icrn = "a + b -> c + d; d + c -> e + f"
        fcrn, fs = parse_crn(fcrn)
        icrn, _ = parse_crn(icrn)

        i01 = {'a': ['B'], 'b': ['A'], 'c': ['C'], 'd': ['D'], 'e': ['C'], 'f': ['D']}
        i02 = {'a': ['B'], 'b': ['A'], 'c': ['C'], 'd': ['D'], 'e': ['D'], 'f': ['C']}
        i03 = {'a': ['B'], 'b': ['A'], 'c': ['C'], 'd': ['D'], 'e': [], 'f': ['C', 'D']}
        i04 = {'a': ['B'], 'b': ['A'], 'c': ['C'], 'd': ['D'], 'e': ['C', 'D'], 'f': []}
        i05 = {'a': ['B'], 'b': ['A'], 'c': ['D'], 'd': ['C'], 'e': ['C'], 'f': ['D']}
        i06 = {'a': ['B'], 'b': ['A'], 'c': ['D'], 'd': ['C'], 'e': ['D'], 'f': ['C']}
        i07 = {'a': ['B'], 'b': ['A'], 'c': ['D'], 'd': ['C'], 'e': ['C', 'D'], 'f': []}
        i08 = {'a': ['B'], 'b': ['A'], 'c': ['D'], 'd': ['C'], 'e': [], 'f': ['C', 'D']}
        i09 = {'a': ['A'], 'b': ['B'], 'c': ['C'], 'd': ['D'], 'e': ['C'], 'f': ['D']}
        i10 = {'a': ['A'], 'b': ['B'], 'c': ['C'], 'd': ['D'], 'e': ['D'], 'f': ['C']}
        i11 = {'a': ['A'], 'b': ['B'], 'c': ['C'], 'd': ['D'], 'e': [], 'f': ['C', 'D']}
        i12 = {'a': ['A'], 'b': ['B'], 'c': ['C'], 'd': ['D'], 'e': ['C', 'D'], 'f': []}
        i13 = {'a': ['A'], 'b': ['B'], 'c': ['D'], 'd': ['C'], 'e': ['C'], 'f': ['D']}
        i14 = {'a': ['A'], 'b': ['B'], 'c': ['D'], 'd': ['C'], 'e': ['D'], 'f': ['C']}
        i15 = {'a': ['A'], 'b': ['B'], 'c': ['D'], 'd': ['C'], 'e': ['C', 'D'], 'f': []}
        i16 = {'a': ['A'], 'b': ['B'], 'c': ['D'], 'd': ['C'], 'e': [], 'f': ['C', 'D']}
        i17 = {'a': ['B'], 'b': ['A'], 'c': ['C', 'D'], 'd': [], 'e': ['C'], 'f': ['D']}
        i18 = {'a': ['B'], 'b': ['A'], 'c': ['C', 'D'], 'd': [], 'e': ['D'], 'f': ['C']}
        i19 = {'a': ['B'], 'b': ['A'], 'c': [], 'd': ['C', 'D'], 'e': ['C'], 'f': ['D']}
        i20 = {'a': ['B'], 'b': ['A'], 'c': [], 'd': ['C', 'D'], 'e': ['D'], 'f': ['C']}
        i21 = {'a': ['A'], 'b': ['B'], 'c': ['C', 'D'], 'd': [], 'e': ['C'], 'f': ['D']}
        i22 = {'a': ['A'], 'b': ['B'], 'c': ['C', 'D'], 'd': [], 'e': ['D'], 'f': ['C']}
        i23 = {'a': ['A'], 'b': ['B'], 'c': [], 'd': ['C', 'D'], 'e': ['C'], 'f': ['D']}
        i24 = {'a': ['A'], 'b': ['B'], 'c': [], 'd': ['C', 'D'], 'e': ['D'], 'f': ['C']}
        bisims = list(crn_bisimulations(fcrn, icrn))
        if len(bisims) != 24:
            print()
            for e, b in enumerate(bisims, 1):
                print(e, b, b in [i01, i02, i03, i04, i05, i06, i07, i08, i09, i10, 
                                  i11, i12, i13, i14, i15, i16, i17, i18, i19, i20, i21, i22, i23, i24])
        assert len(bisims) == 24
        assert i01 in bisims
        assert i02 in bisims
        assert i03 in bisims
        assert i04 in bisims
        assert i05 in bisims
        assert i06 in bisims
        assert i07 in bisims
        assert i08 in bisims
        assert i09 in bisims
        assert i10 in bisims
        assert i11 in bisims
        assert i12 in bisims
        assert i13 in bisims
        assert i14 in bisims
        assert i15 in bisims
        assert i16 in bisims
        assert i17 in bisims
        assert i18 in bisims
        assert i19 in bisims
        assert i20 in bisims
        assert i21 in bisims
        assert i22 in bisims
        assert i23 in bisims
        assert i24 in bisims

    def test_1f_3i(self):
        fcrn = "A + B -> C + D"
        icrn = "a + b <=> c + d; d + c -> e + f"
        fcrn, fs = parse_crn(fcrn)
        icrn, _ = parse_crn(icrn)

        # does not pass permissive.
        bisims = list(crn_bisimulations(fcrn, icrn))
        assert len(bisims) == 0

    def test_2f_2i(self):
        fcrn = "A + B -> C + D; C + D -> E + F"
        icrn = "a + b -> c + d; d + c -> e + f"
        fcrn, fs = parse_crn(fcrn)
        icrn, _ = parse_crn(icrn)

        i01 = {'a': ['B'], 'b': ['A'], 'c': ['C'], 'd': ['D'], 'e': ['F'], 'f': ['E']}
        i02 = {'a': ['B'], 'b': ['A'], 'c': ['C'], 'd': ['D'], 'e': ['E'], 'f': ['F']}
        i03 = {'a': ['B'], 'b': ['A'], 'c': ['D'], 'd': ['C'], 'e': ['F'], 'f': ['E']}
        i04 = {'a': ['B'], 'b': ['A'], 'c': ['D'], 'd': ['C'], 'e': ['E'], 'f': ['F']}
        i05 = {'a': ['A'], 'b': ['B'], 'c': ['C'], 'd': ['D'], 'e': ['F'], 'f': ['E']}
        i06 = {'a': ['A'], 'b': ['B'], 'c': ['C'], 'd': ['D'], 'e': ['E'], 'f': ['F']}
        i07 = {'a': ['A'], 'b': ['B'], 'c': ['D'], 'd': ['C'], 'e': ['F'], 'f': ['E']}
        i08 = {'a': ['A'], 'b': ['B'], 'c': ['D'], 'd': ['C'], 'e': ['E'], 'f': ['F']}
        i09 = {'a': ['B'], 'b': ['A'], 'c': ['D'], 'd': ['C'], 'e': ['F'], 'f': ['E']}
        i10 = {'a': ['A'], 'b': ['B'], 'c': ['D'], 'd': ['C'], 'e': ['F'], 'f': ['E']}
        i11 = {'a': ['B'], 'b': ['A'], 'c': ['D'], 'd': ['C'], 'e': ['E'], 'f': ['F']}
        i12 = {'a': ['A'], 'b': ['B'], 'c': ['D'], 'd': ['C'], 'e': ['E'], 'f': ['F']}
        i13 = {'a': ['B'], 'b': ['A'], 'c': ['C'], 'd': ['D'], 'e': ['F'], 'f': ['E']}
        i14 = {'a': ['A'], 'b': ['B'], 'c': ['C'], 'd': ['D'], 'e': ['F'], 'f': ['E']}
        i15 = {'a': ['B'], 'b': ['A'], 'c': ['C'], 'd': ['D'], 'e': ['E'], 'f': ['F']}
        i16 = {'a': ['A'], 'b': ['B'], 'c': ['C'], 'd': ['D'], 'e': ['E'], 'f': ['F']}

        bisims = list(crn_bisimulations(fcrn, icrn))
        #for e, b in enumerate(bisims, 1):
        #    print(e, {k:v for k, v in sorted(b.items())})
        assert len(bisims) == 16

    def test_order_formals_bug(self):
        fcrn = "A + B -> C"
        icrn = "B_1_ + i7 -> i684 + i17; A <=> i7; i17 -> C_1_ + i29"

        fcrn, fs = parse_crn(fcrn)
        icrn, _ = parse_crn(icrn)

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
        bisims = list(crn_bisimulations(fcrn, icrn))
        #for e, r in enumerate(bisims):
        #    print(f'{e:2d}: {sorted(i.items())}')
        assert len(bisims) == 10
        assert i01 in bisims
        assert i02 in bisims
        assert i03 in bisims
        assert i04 in bisims
        assert i05 in bisims
        assert i06 in bisims
        assert i07 in bisims
        assert i08 in bisims
        assert i09 in bisims
        assert i10 in bisims

@unittest.skipIf(SKIP_DEBUG, "skipping tests for debugging")
class FastBisimulationTests(unittest.TestCase):
    def test_example_01(self):
        # A sample test to aggree on a new interface for bisimulation.
        fcrn = "A->B"
        ecrn = "A<=>i19; i19<=>i39+X; i39->i71+i72"

        fcrn, fs = parse_crn(fcrn)
        ecrn, _ = parse_crn(ecrn)
        partial = {sp: [sp] for sp in fs}

        v, i = crn_bisimulation_test(fcrn, ecrn, fs, interpretation = partial, permissive = 'graphsearch')
        self.assertTrue(v)

        v, i = crn_bisimulation_test(fcrn, ecrn, fs, interpretation = partial, permissive = 'loopsearch')
        self.assertTrue(v)

        v, i = crn_bisimulation_test(fcrn, ecrn, fs, interpretation = partial, permissive = 'bruteforce')
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
                                      permissive = 'bruteforce')
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
                                      permissive = 'bruteforce')
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
                                      permissive = 'bruteforce')
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
        partial = {'x2': ['B', 'D']}
        v, _ = crn_bisimulation_test(fcrn, icrn, fs, interpretation = partial, permissive = 'graphsearch')
        self.assertTrue(v)
        v, _ = crn_bisimulation_test(fcrn, icrn, fs, interpretation = partial, permissive = 'loopsearch')
        self.assertTrue(v)
        v, _ = crn_bisimulation_test(fcrn, icrn, fs, interpretation = partial, permissive = 'bruteforce')
        self.assertTrue(v)

        partial['x3'] = ['C']
        v, _ = crn_bisimulation_test(fcrn, icrn, fs, interpretation = partial, permissive = 'graphsearch')
        self.assertFalse(v)
        v, _ = crn_bisimulation_test(fcrn, icrn, fs, interpretation = partial, permissive = 'loopsearch')
        self.assertFalse(v)
        v, _ = crn_bisimulation_test(fcrn, icrn, fs, interpretation = partial, permissive = 'bruteforce')
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

        v, i1 = crn_bisimulation_test(fcrn, icrn, fs, interpretation = ipart1)
        self.assertTrue(v)
        self.assertDictEqual(i1, ifull1)

        ifull2 = {'B': ['B'],
                  'x1': ['B'],
                  'x2': ['B'],
                  'x3': [],
                  'x4': []}

        ipart2 = {'B': ['B'],
                  'x2': ['B']}

        v, i2 = crn_bisimulation_test(fcrn, icrn, fs, interpretation = ipart2)
        self.assertTrue(v)
        self.assertDictEqual(i2, ifull2)

    def test_example_05(self):
        # Issue fixed: Naming species in certain ways broke bisimulation
        fcrn = "A + C -> A + B"
        icrn = """A <=> x1 + e45 
                C + x1 <=> x3 + x4 
                x3 -> A + B + x5"""

        fcrn, fs = parse_crn(fcrn)
        icrn, _ = parse_crn(icrn)

        inter = {'A': ['A'],
                 'B': ['B'],
                 'C': ['C'],
                 'x1': ['A'],
                 'x3': ['A', 'C']}

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
        assert not v 

    def test_QingDong_crn6_i02_gs_bf(self):
        fcrn, fs = parse_crn('tests/crns/crn6.crn', is_file = True)
        icrn, _ = parse_crn('tests/crns/icrns/crn6_qingdong_thesis.crn', is_file = True)
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
                                     interpretation = inter_02, 
                                     permissive = 'graphsearch')
        self.assertTrue(v)
        v, _ = crn_bisimulation_test(fcrn, icrn, fs, 
                                     interpretation = inter_02, 
                                     permissive = 'bruteforce')
        self.assertTrue(v)

@unittest.skipIf(SKIP_DEBUG, "skipping tests for debugging")
class ModularBisimulationTests(unittest.TestCase):
    def test_qian_roessler_modular(self):
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
        backup = {sp: [sp] for sp in fs}
        v, i = modular_crn_bisimulation_test(fcrns, icrns, fs, partial)
        self.assertTrue(v)
        self.assertDictEqual(partial, backup)

        fcrn, _ = parse_crn('tests/crns/roessler_01.crn', is_file = True)
        icrn, _ = parse_crn('tests/crns/icrns/roessler_qian2011.crn', is_file = True)
        v, i = crn_bisimulation_test(fcrn, icrn, fs, interpretation = i)
        self.assertTrue(v)

    def test_modularity_example_01(self):
        module = """ a <=> i1
                     b + i1 -> i2 + w3
                     i2 -> c + w4
                 """
        module, _ = parse_crn(module)
        fsc = {'A', 'B', 'C'}
        isc = {'a', 'b', 'c'}

        bisim = {'a': ['A'], 'b': ['B'], 'c': ['C'], 'i1': ['A'], 'i2': ['C'], 'w3': [], 'w4': []}
        assert passes_modularity_condition(bisim, module, isc, fsc) is True
        bisim = {'a': ['B'], 'b': ['A'], 'c': ['C'], 'i1': ['B'], 'i2': ['C'], 'w3': [], 'w4': []}
        assert passes_modularity_condition(bisim, module, isc, fsc) is True
        bisim = {'a': ['A'], 'b': ['B'], 'c': ['C'], 'i1': ['A'], 'i2': ['A','B'], 'w3': [], 'w4': []}
        assert passes_modularity_condition(bisim, module, isc, fsc) is False

    def test_modularity_example_02(self):
        module = """ b <=> e1
                     e1 -> e2 + e3 + e4
                     e4 -> e1 + e5
                 """
        module, _ = parse_crn(module)
        fsc = {'B'}
        isc = {'b'}

        minter = {'b': ['B'], 'e1': ['B'], 'e2': [], 'e3': [], 'e4': [], 'e5': []} 
        assert passes_modularity_condition(minter, module, isc, fsc) is True

        minter = {'b': ['B'], 'e1': ['B'], 'e2': [], 'e3': [], 'e4': ['B'], 'e5': []} 
        assert passes_modularity_condition(minter, module, isc, fsc) is True

        minter = {'b': ['B'], 'e1': ['B', 'B'], 'e2': [], 'e3': [], 'e4': [], 'e5': []} 
        assert passes_modularity_condition(minter, module, isc, fsc) is False

        minter = {'b': ['B'], 'e1': ['B'], 'e2': ['B'], 'e3': [], 'e4': [], 'e5': []} 
        assert passes_modularity_condition(minter, module, isc, fsc) is False
        isc = {'b', 'e2'}
        assert passes_modularity_condition(minter, module, isc, fsc) is True
        isc = {'b', 'e5'}
        minter = {'b': ['B'], 'e1': [], 'e2': [], 'e3': [], 'e4': [], 'e5': ['B']} 
        assert passes_modularity_condition(minter, module, isc, fsc) is True

        minter = {'b': ['B'], 'e1': ['B'], 'e2': ['B'], 'e3': [], 'e4': [], 'e5': ['B']} 
        assert passes_modularity_condition(minter, module, isc, fsc) is False

    def test_modularity_example_03(self):
        module = """
                    a <=> e6
                    a <=> e1
                    a + e6 <=> e1
                """
        module, _ = parse_crn(module)
        fsc, isc = {'A'}, {'a'}
        inter = {'a': ['A'], 'e1': ['A', 'A'], 'e6': ['A']}
        assert passes_modularity_condition(inter, module, isc, fsc) is True

@unittest.skipIf(SKIP_SLOW or SKIP_DEBUG, "skipping tests for debugging")
class SlowBisimulationTests(unittest.TestCase):
    def test_QingDong_crn6_i1_gs(self):
        # NOTE: around 3 minutes
        fcrn, fs = parse_crn('tests/crns/crn6.crn', is_file = True)
        icrn, _ = parse_crn('tests/crns/icrns/crn6_qingdong_thesis.crn', is_file = True)
        inter_01 = {'i778': ['Y'],
                    'i575': ['X'],
                    'i599': ['C'],
                    'i2232': ['A'],
                    'i73': ['B']}

        v, _ = crn_bisimulation_test(fcrn, icrn, fs, 
                                     interpretation = inter_01,
                                     permissive = 'graphsearch')
        self.assertTrue(v)

    def test_QingDong_crn6_i1_bf(self):
        # NOTE: around 3 minutes
        fcrn, fs = parse_crn('tests/crns/crn6.crn', is_file = True)
        icrn, _ = parse_crn('tests/crns/icrns/crn6_qingdong_thesis.crn', is_file = True)
        inter_01 = {'i778': ['Y'],
                    'i575': ['X'],
                    'i599': ['C'],
                    'i2232': ['A'],
                    'i73': ['B']}
        v, _ = crn_bisimulation_test(fcrn, icrn, fs, 
                                     interpretation = inter_01, 
                                     permissive = 'bruteforce')
        self.assertTrue(v)

    def test_QingDong_crn6_gs(self):
        # NOTE: less than 5 min.
        fcrn, fs = parse_crn('tests/crns/crn6.crn', is_file = True)
        icrn, _ = parse_crn('tests/crns/icrns/crn6_qingdong_thesis.crn', is_file = True)
        # These tests complete in less than 10 minutes
        v, _ = crn_bisimulation_test(fcrn, icrn, fs, permissive = 'graphsearch')
        self.assertTrue(v)

    def test_QingDong_crn6_bf(self):
        # NOTE: less than 10 min.
        fcrn, fs = parse_crn('tests/crns/crn6.crn', is_file = True)
        icrn, _ = parse_crn('tests/crns/icrns/crn6_qingdong_thesis.crn', is_file = True)
        v, _ = crn_bisimulation_test(fcrn, icrn, fs, permissive = 'bruteforce')
        self.assertTrue(v)

    def dtest_QingDong_crn6_i02_ls(self):
        # TODO: does not finish ... probably a bug!
        fcrn, fs = parse_crn('tests/crns/crn6.crn', is_file = True)
        icrn, _ = parse_crn('tests/crns/icrns/crn6_qingdong_thesis.crn', is_file = True)
        inter_02 = {'i842': ['Y', 'X', 'A'],
                    'i394': ['X', 'Y', 'X'],
                    'i119': ['X', 'B', 'A'],
                    'i2300': ['A', 'C'],
                    'i778': ['Y'],
                    'i575': ['X'],
                    'i599': ['C'],
                    'i2232': ['A'],
                    'i73': ['B']}
        v, _ = crn_bisimulation_test(fcrn, icrn, fs, interpretation = inter_02, permissive = 'loopsearch')
        self.assertTrue(v)

    def dtest_QingDong_crn6_i1_ls(self):
        # TODO: does not finish ... probably a bug!
        fcrn, fs = parse_crn('tests/crns/crn6.crn', is_file = True)
        icrn, _ = parse_crn('tests/crns/icrns/crn6_qingdong_thesis.crn', is_file = True)
        inter_01 = {'i778': ['Y'],
                    'i575': ['X'],
                    'i599': ['C'],
                    'i2232': ['A'],
                    'i73': ['B']}
        v, _ = crn_bisimulation_test(fcrn, icrn, fs, 
                                     interpretation = inter_01, 
                                     permissive = 'loopsearch')
        self.assertTrue(v)

    def dtest_QingDong_crn6_ls(self):
        # TODO: test again after i1 and i2 terminate.
        fcrn, fs = parse_crn('tests/crns/crn6.crn', is_file = True)
        icrn, _ = parse_crn('tests/crns/icrns/crn6_qingdong_thesis.crn', is_file = True)
        v, _ = crn_bisimulation_test(fcrn, icrn, fs, permissive = 'loopsearch')
        self.assertTrue(v)

    def dtest_qian_roessler_full(self):
        # TODO: KILLED / OUT OF MEMORY
        # Seems like "space-efficient still takes quite some memory, so we'll
        # have to debug there.
        (fcrn, fs) = parse_crn('tests/crns/roessler_01.crn', is_file = True)
        (icrn, _) = parse_crn('tests/crns/icrns/roessler_qian2011.crn', is_file = True)
        partial = {sp: [sp] for sp in fs}
        v, i = crn_bisimulation_test(fcrn, icrn, fs, partial, searchmode = 'space-efficient')
        self.assertTrue(v)


if __name__ == '__main__':
    unittest.main()
