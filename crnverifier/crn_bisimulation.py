#
#  crnverifier/crn_bisimulation.py
#  Original source from the Nuskell compiler project
#
#  Authors:
#   - Qing Dong
#   - Robert Johnson
#   - Stefan Badelt
#
import logging
log = logging.getLogger(__name__)

import copy
import math
from collections import Counter
from itertools import product, permutations, combinations, chain

from .utils import pretty_rxn
from .utils import interpret as interpretL
from .deprecated import subsets, enum 

class SpeciesAssignmentError(Exception):
    pass

class SearchDepthExceeded(Exception):
    pass

class CRNBisimulationError(Exception):
    pass

# Conversion between list and counter objects.
def rl(rxn):
    return [list(part.elements()) for part in rxn]

def rc(rxn):
    return [Counter(part) for part in rxn]

def inter_counter(inter):
    # For internal representation of interpretations 
    # return interpretations of the format: inter[str] = Counter
    return {k: Counter(v) for k, v in inter.items()}

def inter_list(inter):
    # For consistency with the rest of the package, 
    # return interpretations of the format: inter[str] = list
    return {k: list(v.elements()) for k, v in inter.items()}

# Conversion between internal and external representation of implementation species.
def deformalize(k, fs):
    return f'i{{{k}}}' if k in fs else k

def formalize(intrp):
    intr = {}
    for sp in intrp:
        if sp[:2] == 'i{':
            isp = sp[2:-1]
            intr[isp] = intrp[sp]
        else:
            intr[sp] = intrp[sp]
    return intr

# Utils for list based multiset operations.
def subsetsL(x):
    """ Generate all (uniqe) multi-subsets of a list x. """
    return chain(*[combinations(x, l) for l in range(len(x), -1, -1)])

def subtractL(l1, l2, strict = True):
    """ Returns a new list: l1 - l2. """
    l1 = l1[:]
    for el in l2:
        try:
            l1.remove(el)
        except ValueError as err:
            if strict:
                raise err
    return l1

def is_contained(a, b):
    """ bool: True if list a is fully contained in list b. """
    b = b[:]
    try:
        [b.remove(s) for s in a]
    except ValueError as err:
        return False
    return True

def enumL(n, l, weights = None):
    """ Returns combinations to assign elements of list l to n variables. """
    if n == 0:
        yield []
        return

    if weights is None:
        weights = [1 for x in range(n)]
    w = weights[0]

    def ldiv(l, s):
        """ True, if all elements are multiples of s """
        if s == 1 or len(l) == 0:
            return l
        l = Counter(l)
        for k, v in l.items():
            if v % s:
                break
            l[k] = int(v/s)
        else:
            return tuple(l.elements())
        return None
        
    if n == 1:
        l = ldiv(l, w)
        if l is None:
            raise SpeciesAssignmentError
        yield [tuple(l)]
        return

    for ss in set(subsetsL(l)):
        assert is_contained(ss, l)
        wss = ldiv(ss, w)
        if wss is None:
            continue
        try:
            for j in enumL(n-1, subtractL(l, ss), weights[1:]):
                yield [wss] + j
        except SpeciesAssignmentError:
            continue
    return

def msleq(x,y):
    # True if (multisets) x <= y (vector comparison)
    for k in x:
        if x[k] > y[k]:
            return False
    return True

def interpret(s, intrp):
    # return interpretation of s according to intrp
    ss = s.copy()
    for k in list(s):
        if k in intrp:
            v = ss.pop(k)
            for i in range(v):
                ss += intrp[k]
    return ss

def interleq(x, y, intrp):
    # True if m(x) <= m(y) with m given by interpretation intrp
    return msleq(interpret(x,intrp),interpret(y,intrp))

def subst(crn, intrp, counter = True):
    """ Substitute implementation species with formal species in CRN.  """
    if counter:
        return [[interpret(j, intrp) for j in rxn] for rxn in crn]
    else:
        return [[interpretL(part, intrp) for part in rxn] for rxn in crn]

def same_reaction(irxn, frxn, fs, counter = True):
    """ Check if irxn *could* implement frxn.

    This assumes that irxn is already interpreted using the formal species fs.
    SB: Note, that I introduced a new expression expression here, see below.
    """
    if counter:
        fr = set(frxn[0] - irxn[0]) # exess formal reactants 
        ir = set(irxn[0] - frxn[0]) # exess implementation reactants 
        fp = set(frxn[1] - irxn[1]) # exess formal products
        ip = set(irxn[1] - frxn[1]) # exess implementation products
    else:
        fr = set(subtractL(frxn[0], irxn[0], strict = False)) # exess formal reactants 
        ir = set(subtractL(irxn[0], frxn[0], strict = False)) # exess implementation reactants 
        fp = set(subtractL(frxn[1], irxn[1], strict = False)) # exess formal products
        ip = set(subtractL(irxn[1], frxn[1], strict = False)) # exess implementation products
    if ir & fs or ip & fs:
        # There are excess formal reactants or excess formal products
        # in the implementation reaction, so this reaction cannot
        # implement the formal reaction.
        return False
    if len(fr) and len(ir) == 0:
        # There are excess formal reactants and 
        # no more implementation reactants.
        return False 
    if len(fp) and len(ip) == 0:
        # There are excess formal products and 
        # no more implementation products.
        return False
    if len(fp) and len(ir) and len(ip) and ir == ip:
        # NOTE: This line is new:
        # Example:    A -> B + C
        # cannot be:  A + y -> B + y
        return False
    return True

def updateT(fcrn, icrn, fs, counter = True):
    """ Calculate a table with bool entries.

    Assumes an implementation CRN where all implementation species
    have been replaced with their formal interpretation. Thus, the
    implementation CRN contains a mix of formal species and/or 
    implementation species.

    row num = len(icrn)
    col num = len(fcrn) + 1

    """
    # RJ: This should be logically equivalent to the UpdateTable in the MS
    # thesis for compiled DNA reactions (we think).
    # WARNING: If an implementation CRN has directly catalytic species, the
    # code below may fail (though thesis psuedocode is correct).  
    # E.g.  i3 + i7 --> i12 + i7 + i8
    r = []
    for irxn in icrn:
        rr = [same_reaction(irxn, frxn, fs, counter) for frxn in fcrn]

        if counter:
            ir = set(irxn[0] - irxn[1]) # "non-trivial" reactants
            ip = set(irxn[1] - irxn[0]) # "non-trivial" products
        else:
            ir = set(subtractL(irxn[0], irxn[1], strict = False)) # "non-trivial" reactants
            ip = set(subtractL(irxn[1], irxn[0], strict = False)) # "non-trivial" products
        if (ir & fs and ip <= fs) or (ip & fs and ir <= fs):
            # If the implementation reactants contain a formal species and
            # all implementation products are (different) formal species, or
            # if the implementation products contain a formal species and
            # all implementation reactants are (different) formal species:
            rr.append(False)
        else: # Could be trival
            rr.append(True)
        r.append(rr)
    return r

def checkT(T):
    """ Check (partial) interpretation for the delimiting condition.

    If there is a whole row, or a whole (non-trivial) column False, then 
    the permissive condition is violated and we have to roll back.

    Returns:
        bool: True if there is no row or (non-trivial) colum with all False values.
    """
    for row in T:
        if all(not b for b in row):
            return False
    n = len(T[0])
    for e, col in enumerate(zip(*T), 1):
        if e != n and all(not b for b in col):
            return False
    return True

def formal_states(fstate, inter, counter = False):
    """ Generate a minimal set of istates which interpret to fstate.

    The returned implementation states may interpret to supersets of fstate.
    """
    if len(fstate) == 0:
        yield [] if not counter else Counter()
    else:
        fs = set(fstate).pop()



        fs1 = fstate[0]
        for k in inter:
            if fs1 in inter[k]:
                for out in formal_states(subtractL(fstate, inter[k], strict = False), 
                                                    inter,
                                                    counter):
                    yield ([k] + out) if not counter else (Counter({k:1}) + out)

def passes_modularity_condition(bisim, icrn, common_is, common_fs):
    """ Check if a bisimulation satisfies the modularity condition.

    The modularity condition is formulated with respect to a subset of *common
    implementation* and *common formal* species, where each common
    implementation species must interpet to a common formal species.  The
    bisimulation is modular if every implementation species can turn into
    common implementation species with the same interpretation via trivial
    reactions (that is the interpretation of the reaction must be trival).
    For details check JDW2019: Definition 3.3 (Modularity condition).

    SB: contains a slight change in iteration compared to "theOne" before.

    (This is basically the graphsearch algorithm to check the permissive
    condition, where all "minimal states" are each exactly one implementation
    species.)

    Args:
        bisim (dict): An interpretation which is a CRN bisimulation.
        icrn (list[lists]): An implementation CRN.
        common_is (set): A set of common implementation species.
        common_fs (set): A set of common formal species.

    Returns:
        bool: True if modularity condition is satisfied.
    """
    # All common implementation species are part of the interpretation.
    assert all(ci in bisim for ci in common_is)
    # All interpretations of common implementation species 
    # are part of common formal species.
    # assert all(cf in common_fs for ci in common_is for cf in bisim[ci])

    def Y(s):
        return s in common_is
    def Z(s):
        return len(set(bisim[s]) & common_fs) == 0

    trivial_rxns = [] # trivial reactions
    for rxn in icrn:
        iR = interpretL(rxn[0], bisim) 
        iP = interpretL(rxn[1], bisim) 
        if sorted(iR) == sorted(iP):
            trivial_rxns.append(rxn)

    done = {s for s in bisim if Y(s) or Z(s)}
    todo = {s: [set(), set()] for s in bisim if s not in done} 
    # todo[s] = [reach, produce]

    reset = len(todo) > 0
    while reset:
        reset = False
        for r in list(todo.keys()): # start at any of the non-common intermediates
            [r_reach, r_produce] = todo[r]
            for R, P in trivial_rxns:
                if r not in R:
                    continue
                nR = R[:]
                nR.remove(r)
                if set(nR) > r_produce:
                    # There is at least one more reactant which (as far as we
                    # have seen) cannot be produced by r.
                    continue

                nullsp, mystic = set(), set()
                for p in P:
                    if p in done or (not is_contained(bisim[r], bisim[p])): 
                        # If p is done, then any other p in P is either {} or
                        # it is also a reactant which must be checked
                        # separately. Alternatively, if the interpretation
                        # of products contains the interpretation
                        # of reactants and more, then we can call this one
                        # done and move on to check p.
                        done.add(r)
                        del todo[r]
                        break
                    elif len(bisim[p]) == 0:
                        nullsp.add(p)
                    else:
                        mystic.add(p)
                if r in done:
                    reset = True
                    break
                # Now investigate the mysterious "non-common implementation"
                # products which interpret to a formal species, but are not
                # known to be exit states.
                for p in mystic: 
                    if p not in r_reach:
                        r_reach.add(p)
                        reset = True
                    [p_reach, p_produce] = todo[p]
                    if not (p_reach <= r_reach):
                        # anything reachable by products is reachable by reactants.
                        r_reach |= p_reach
                        reset = True
                    if r in p_reach:
                        # we know r -> p -> r, so all r_nullsp and p_nullsp 
                        # can be produced at infinite amounts.
                        loopable = nullsp | p_produce
                        if not (loopable <= r_produce):
                            r_produce |= loopable
                            reset = True
    return len(todo) == 0

def check_permissive_reactionsearch(fcrn, icrn, fs, intrp):
    """ Check the permissive condition.
    
    Uses the original formal CRN  and original implementation CRN
    (without substitutions).

    Args:
        fcrn: The original formal CRN
        icrn: The original implementation CRN
        fs: The formal species.

    """
    permissive_depth = None
    sicrn = subst(icrn, intrp)
    T = updateT(fcrn, sicrn, fs)
    assert checkT(T) # Assuming this has been checked before calling permissive.

    nulls = [k for k in intrp if not len(list(intrp[k]))]
    trivial_rxns = [irxn for e, irxn in enumerate(icrn) if T[e][-1]]

    def cnstr(s):
        """
        # given formal state s, generate minimal set of impl. states which
        #  interpret to a superset of s
        # = concatenation of any x such that s[0] in m(x) with
        #  all minimal implementations of s - m(x)
        """
        if len(s) == 0:
            yield Counter()
        else:
            s1 = list(s.keys())[0]
            for k in intrp:
                if s1 in intrp[k]:
                    assert not (s - intrp[k])[s1] >= s[s1]
                    for out in cnstr(s - intrp[k]):
                        yield Counter({k:1}) + out

    def reactionsearch(s, d):
        # s is the implementation state, d is the depth.
        # try to find (via trivial reactions) a state in which a reaction in fr can fire.
        # fr is a particular formal reaction along with all implementation reactions that interpret to it. 
        # tr is the list of all trivial reactions in the implementation.
        # fail if depth d of trivial reaction steps is exceeded without firing a reaction in fr.
        if permissive_depth and d > permissive_depth:
            return None
        hashee = list(s.elements())
        for i in range(len(hashee)):
            try:
                if hashee[i][0] == 'impl':
                    hashee[i] = hashee[i][1]
            except (IndexError, TypeError):
                pass
        hashee = tuple(sorted(hashee))
        if hashee in hasht:
            return False
        else:
            hasht.add(hashee)
        for i in fr:
            if list((i[0] - s).keys()) == []:
                return True
        ret = False
        for i in trivial_rxns:
            if list((i[0] - s).keys()) == []:
                t = (s - i[0]) + i[1]
                out = reactionsearch(t, d+1)
                if out:
                    return True
                elif out is None:
                    ret = None
        return ret

    for i, frxn in enumerate(fcrn):
        # build fr for this formal reaction
        fr = [irxn for e, irxn in enumerate(icrn) if T[e][i]]

        ist = cnstr(frxn[0])
        # At this point, "ist" contains an exhaustive and sufficient list of
        # possible initial implementation states in which the current formal
        # reaction #i must be able to fire (via a series of trivial reactions),
        # Note that we will only want to test states in "ist" that interpret to
        # a state in which #i can fire.

        tested = []  # to avoid testing a superset of some state that's already been tested
        spaceefficient = True # ... but only if we have space to store them
        for j in ist:
            tmp = interpret(j,intrp)
            if msleq(frxn[0], tmp):  # only test if reactants j interpret to allow #i to fire
                t = False
                for k in tested:
                    # will be a 0-length loop if spaceefficient
                    if msleq(k, j):
                        t = True
                        break
                if t:
                    continue
                hasht = set([])
                found = reactionsearch(j, 0)
                if not found: 
                    return False, (frxn, j)
                elif not spaceefficient:
                    tested.append(j)
    return True, 0

def check_permissive_loopsearch(fcrn, icrn, fs, inter):
    """ The 'loopsearch' algorithm to check the permissive condition. 
    """
    sicrn = subst(icrn, inter, counter = False)
    T = updateT(fcrn, sicrn, fs, counter = False)
    assert checkT(T) # Assuming this has been checked before calling permissive.
    assert all(fs | set(v) == fs for k, v in inter.items())

    log.debug(f'The implementation CRN:\n' + '\n'.join(
        [f'{e} {t=} {pretty_rxn(sicrn[e])} ({pretty_rxn(icrn[e])})' \
                for e, t in enumerate(T)]))
    log.debug(f'{inter=}')

    # Potential null species.
    nulls = [k for k, v in inter.items() if len(v) == 0]

    # Those are all remaining trivial reactions.
    trivial_rxns = [irxn for e, irxn in enumerate(icrn) if T[e][-1]]

    def midsearch(start, goal, pickup, ignore, fstates, k):
        # search for a path from start to goal (states of non-null species)
        #  which produces at least pickup null species
        #  assuming it already has infinite copies of everything in ignore
        #  of length at most 2^k
        # if goal is None, the goal is any reaction in fr
        if goal is not None:
            if ignore.issuperset(subtractL(list(goal), start, False)):
                return True
            if not is_contained(interpretL(goal, inter), 
                                interpretL(start, inter)):
                return False

        if k == 0:
            if goal is None:
                for rx in fr:
                    if ignore.issuperset(subtractL(rx[0], start, False)):
                        return True
                return False
            
            for rx in trivial_rxns:
                if ignore.issuperset(subtractL(rx[0], start, False)):
                    # every element of rx[0] - start (multiset)
                    #  is also an element of ignore (set)
                    # e.g. true on rx[0] - start = {|a,a,b|}, ignore = {a,b,c}
                    if is_contained(goal, subtractL(start, rx[0], False) + rx[1]):
                        return True

        else:
            if midsearch(start, goal, pickup, ignore, fstates, k-1):
                return True
            for part in subsetsL(pickup):
                for mid in fstates:
                    if midsearch(start, mid, set(part), ignore, fstates, k-1) \
                       and ((not is_contained(interpretL(start, inter), interpretL(mid, inter))) or
                            midsearch(mid,goal,pickup-set(part),ignore, fstates,k-1)):
                        return True

        return False

    def loopsearch(start, fstates):
        # search for a path from start to any reaction in fr
        rounds, roundequiv = 0, 1
        for nequiv, point in enumerate(fstates, 1):
            if nequiv > roundequiv:
                rounds += 1
                roundequiv *= 2

        for parti in enumL(len(nulls) + 1, nulls):
            log.debug(f'{parti=}')
            part = list(map(set, parti))
            if any([part[i] != set() and part[i+1] == set() for i in range(len(part) - 2)]):
                continue # avoid redundancy

            pickups = [_f for _f in part[:-1] if _f is not None]
            place = start
            ignore = set()
            for pickup in pickups:
                for base in fstates:
                    if midsearch(place, base, set(), ignore, fstates, rounds):
                        if not is_contained(interpretL(place, inter), 
                                            interpretL(base, inter)):
                            return True
                        elif midsearch(base, base, pickup, ignore, fstates, rounds):
                            place = base
                            break
                else:
                    break
                ignore |= pickup
            else:
                if midsearch(place, None, set(), ignore, fstates, rounds):
                    return True
        return False

    mode = 'time-efficient'
    for i, frxn in enumerate(fcrn):
        # build fr for this formal reaction
        fr = [irxn for e, irxn in enumerate(icrn) if T[e][i]]

        fstates = list(map(sorted, formal_states(frxn[0], inter)))
        tested = []  # to avoid testing a superset of some state that's already been tested
        for fstate in fstates:
            # The formal reaction can happen.
            if is_contained(frxn[0], interpretL(fstate, inter)):
                for k in tested:
                    # TODO: might make sense to sort those states, right?
                    if is_contained(k, list(fstate)):
                        break
                else:
                    hasht = set([])
                    if loopsearch(fstate, fstates):
                        if mode == 'time-efficient':
                            tested.append(fstate)
                    else:
                        return False, (frxn, fstate)
    return True, 0

def check_permissive_graphsearch(fcrn, icrn, fs, inter):
    """ The 'graphsearch' algorithm to check the permissive condition. 

    Checks whether every implementation of a formal state that can react in the
    formal CRN by reaction r can also react in the implementation CRN with a 
    reaction that interprets to r.

    """

    sicrn = subst(icrn, inter, counter = False)
    T = updateT(fcrn, sicrn, fs, counter = False)
    assert checkT(T) # assuming this has been checked before calling permissive.
    assert all(fs | set(v) == fs for k, v in inter.items())

    log.debug(f'The implementation CRN:\n' + '\n'.join(
        [f'{e} {t=} {pretty_rxn(sicrn[e])} ({pretty_rxn(icrn[e])})' \
                for e, t in enumerate(T)]))
    log.debug(f'{inter=}')

    # Potential null species.
    nulls = set([k for k, v in inter.items() if len(v) == 0])

    # Those are all remaining trivial reactions.
    trivial_rxns = [irxn for e, irxn in enumerate(icrn) if T[e][-1]]

    max_depth = 0
    for i, frxn in enumerate(fcrn):
        # The implementation reactions for this formal reaction.
        implementations = [irxn for e, irxn in enumerate(icrn) if T[e][i]]
        # All implementation states which must permit the formal reaction.
        fstates = list(map(tuple, map(sorted, formal_states(frxn[0], inter))))
        done = set() # set of species known to implement the current formal rxn.
        todo = {fstate: [set(), set()] for fstate in fstates}
        # todo[fstate] = [set(nulls), set(reachable)]

        log.debug(f'Testing formal reaction "{pretty_rxn(frxn)}" with {fstates=}')
        changed, depth = True, 0
        while changed:
            changed, depth = False, depth + 1
            for fstate in fstates:
                if fstate in done:
                    continue
                for irxn in implementations:
                    # If the implementation state can implement the formal reaction.
                    if set(subtractL(irxn[0], fstate, False)) <= todo[fstate][0]: 
                        done.add(fstate)
                        del todo[fstate]
                        changed = True
                        break
                if fstate in done:
                    continue
                for irxn in trivial_rxns:
                    # If the implementation state can implement the trivial reaction.
                    if set(subtractL(irxn[0], fstate, False)) <= todo[fstate][0]: 
                        # remove the reactants from the implementation state.
                        nstate = subtractL(list(fstate), irxn[0], False) + irxn[1]
                        for fstate2 in fstates:
                            # if one of the other implementation states is reachable
                            if is_contained(fstate2, nstate):
                                if fstate2 in done:
                                    done.add(fstate)
                                    del todo[fstate]
                                    changed = True
                                    break
                                if fstate in todo[fstate2][1]: # I assume its a loop?
                                    s = todo[fstate2][0] | (set(nstate) & nulls)
                                    if not s <= todo[fstate][0]:
                                        todo[fstate][0] |= s
                                        changed = True
                                if fstate2 not in todo[fstate][1]:
                                    todo[fstate][1].add(fstate2)
                                    changed = True
                                if not (todo[fstate2][1] <= todo[fstate][1]):
                                    todo[fstate][1] |= todo[fstate2][1]
                                    changed = True
                    if fstate in done:
                        break
        if todo: 
            return False, fstates
        if max_depth < depth: 
            max_depth = depth
    return True, max_depth
        
def passes_permissive_condition(fcrn, icrn, fs, intrp, permcheck = 'graphsearch'):
    """
    """
    log.debug(f'Checking permissive condition using {permcheck=}.')

    if permcheck == 'graphsearch':
        # Convert to list format.
        fcrn = [rl(frxn) for frxn in fcrn]
        icrn = [rl(irxn) for irxn in icrn]
        inter = inter_list(intrp)
        passes, info = check_permissive_graphsearch(fcrn, icrn, fs, inter)
    elif permcheck == 'loopsearch':
        if False:
            print(f'WARNING: {permcheck=} has a bug (old version).')
            from .deprecated import permissive
            passes, info = permissive(fcrn, icrn, fs, intrp, permcheck)
        else:
            print(f'WARNING: {permcheck=} has a bug.')
            fcrn = [rl(frxn) for frxn in fcrn]
            icrn = [rl(irxn) for irxn in icrn]
            inter = inter_list(intrp)
            passes, info = check_permissive_loopsearch(fcrn, icrn, fs, inter)
    elif permcheck == 'reactionsearch':
        passes, info = check_permissive_reactionsearch(fcrn, icrn, fs, intrp)
    else:
        raise SystemExit(f'Unknown algorithm for permissive condition: {permcheck}')

    return passes, info

def passes_delimiting_condition(fcrn, icrn, fs, inter):
    """ Tests an interpretation for the delimiting condition.

    Note that this function is actually never used, because typically you want
    to keep the table T. In any case, it is tested and yields the correct
    result.
    """
    sicrn = subst(icrn, inter)
    T = updateT(fcrn, sicrn, fs)
    return checkT(T)

def passes_atomic_condition(inter, fs):
    """ Tests an interpretation for the atomic condition.

    For every formal species there exists an implementation species which
    interprets to it.  
    """
    return all(any(set(f) == set(v.elements()) for v in inter.values()) for f in fs)

def search_row(fcrn, icrn, fs, intrp, moves = None, depth = 0, 
               mode = 'time-efficient', permcheck = 'graphsearch'):
    """ Find full interpretations matching every irxn to one frxn or trxn.

    This "row search" finds all valid combinations of 
        implementation reaction -> formal reaction or 
        implementation reaction -> trivial reaction.
    Typically this function is called after search_column, which means we
    already have delimiting and atomic condition satisfied. However, the
    delimiting condition may still break when an additional implementation
    reaction is introduced. 

    """
    assert mode in ('time-efficient', 'space-efficient')
    if moves is None:
        moves = set()

    log.info(f'Searching row at {depth=}'.center(80-(2*depth), '~'))
    log.debug(f'Partial interpretation {inter_list(intrp)=}')
    log.debug(f'Forbidden moves {moves=}')

    sicrn = subst(icrn, intrp)
    T = updateT(fcrn, sicrn, fs)
    if not checkT(T):
        log.debug(f'Delimiting condition not satisfied.')
        raise SpeciesAssignmentError('Delimiting condition not satisfied.')

    # The row indices where unassigned species can be found, but we 
    # remove duplicates of implementation reactions that solve for the
    # same set of implementation species!
    unknown = {tuple(sorted((set(sirxn[0]) | set(sirxn[1])) - fs)): i \
                for i, sirxn in enumerate(sicrn)} 
    unknown = [v for k, v in unknown.items() if len(k)]

    if unknown == []:
        if not passes_atomic_condition(intrp, fs):
            log.debug(f'Atomic condition not satisfied.')
            raise SpeciesAssignmentError('Atomic condition not satisfied.')
        correct, info = passes_permissive_condition(fcrn, icrn, fs, intrp, permcheck)
        if correct:
            log.debug(f'Permissive condition satisfied ({info=}).')
            yield (intrp, moves) if depth > 0 else intrp
            return
        log.debug(f'Permissive condition not satisfied ({info=}).')
        raise SpeciesAssignmentError('Permissive condition not satisfied.')

    unknown = sorted(unknown, key = lambda x: (T[x][-1], T[x][:-1].count(True)))
    # Start with the row with the minimal number of True's (excluding trivial reactions)
    log.debug(f'{unknown=} \n' + '\n'.join(
        [f'{e} {t=} {pretty_rxn(rl(sicrn[e]))} ({pretty_rxn(rl(icrn[e]))})' \
                for e, t in enumerate(T)]))

    later = dict() # less promising partial interpretations ...
    def search(later, moves, n = 0, m = None):
        # use this to define a search hierarchy when sending results into the
        # next round. Currently, n and m are start and stop in a range of 
        # how many species interpet to a single species.
        while later:
            if m and n > m:
                break
            if n in later:
                for (inext, depth, imove) in later[n]:
                    if imove in moves:
                        continue
                    if mode == 'time-efficient':
                        moves.add(imove)
                    try:
                        for ir, mr in search_row(fcrn, icrn, fs, inext, 
                                                 moves = moves, 
                                                 depth = depth + 1, 
                                                 mode = mode, 
                                                 permcheck = permcheck):

                            moves |= mr
                            yield ir, mr
                    except SpeciesAssignmentError:
                        continue
                del later[n]
            n += 1

    for k in unknown:
        irxn = sicrn[k] 

        if T[k][-1] is True: # Assign a trivial reaction.
            log.debug(f'Interpret: {pretty_rxn(rl(irxn))} => trivial')
            fl = Counter({k:v for k, v in (irxn[0] - irxn[1]).items() if k in fs})
            fr = Counter({k:v for k, v in (irxn[1] - irxn[0]).items() if k in fs})
            ul = Counter({k:v for k, v in (irxn[0] - irxn[1]).items() if k not in fs})
            ur = Counter({k:v for k, v in (irxn[1] - irxn[0]).items() if k not in fs})
            if len(ul) or len(ur):
                [kl, vl] = zip(*ul.items()) if len(ul) else [[], []]
                tmpl = enumL(len(ul), list(fr.elements()), vl)
                [kr, vr] = zip(*ur.items()) if len(ur) else [[], []]
                tmpr = enumL(len(ur), list(fl.elements()), vr)

                for i, j in product(tmpl, tmpr):
                    intrpleft = {k: tuple(sorted(v)) for k, v in zip(kl, i)}
                    intrpright = {k: tuple(sorted(v)) for k, v in zip(kr, j)}
                    for key in set(intrpleft) & set(intrpright):
                        if any([intrpleft[key][fsp] != \
                                intrpright[key][fsp] for fsp in fs]):
                            # Incompatible dictionaries!
                            break
                    else:
                        inext = intrp.copy()
                        imove = tuple(sorted(chain(intrpleft.items(), intrpright.items())))
                        inext.update({k: Counter(v) for k, v in imove})
                        level = max(len(v) for k, v in imove)
                        if level in later:
                            later[level].append((inext, depth, imove))
                        else:
                            later[level] = [(inext, depth, imove)]

        for c, frxn in enumerate(fcrn): # Assign a formal reaction.
            if not T[k][c]:
                continue
            log.debug(f'Interpret: {pretty_rxn(rl(irxn))} => {pretty_rxn(rl(frxn))}')
            # left 
            ul = irxn[0] - frxn[0]
            sl = frxn[0] - irxn[0]
            # right
            ur = irxn[1] - frxn[1]
            sr = frxn[1] - irxn[1]

            if len(ul) or len(ur):
                [kl, vl] = zip(*ul.items()) if len(ul) else [[], []]
                [kr, vr] = zip(*ur.items()) if len(ur) else [[], []]
                tmpl = enumL(len(ul), list(sl.elements()), vl)
                tmpr = enumL(len(ur), list(sr.elements()), vr)

                for i, j in product(tmpl, tmpr):
                    intrpleft = {k: tuple(sorted(v)) for k, v in zip(kl, i)}
                    intrpright = {k: tuple(sorted(v)) for k, v in zip(kr, j)}
                    for key in set(intrpleft) & set(intrpright):
                        if any([intrpleft[key][fsp] != \
                                intrpright[key][fsp] for fsp in fs]):
                            # Incompatible dictionaries!
                            break
                    else:
                        inext = intrp.copy()
                        imove = tuple(sorted(chain(intrpleft.items(), intrpright.items())))
                        inext.update({k: Counter(v) for k, v in imove})
                        level = max(len(v) for k, v in imove)
                        if level in later:
                            later[level].append((inext, depth, imove))
                        else:
                            later[level] = [(inext, depth, imove)]
            for i, m in search(later, moves, n = 0, m = 1):
                yield (i, m) if depth > 0 else i
        for i, m in search(later, moves, n = 0, m = 2):
            yield (i, m) if depth > 0 else i
    for i, m in search(later, moves, n = 2):
        yield (i, m) if depth > 0 else i
    assert len(later) == 0
    return

def search_column(fcrn, icrn, fs = None, intrp = None, unknown = None, depth = 0):
    """ Find all interpretation seeds matching every frxn to one irxn.

    This "column search" finds all valid combinations of formal reaction to
    implementation reaction. For example, if there is one formal reaction and
    three compatible implementation reactions, then it will return three
    interpretations that have exactly one reaction interpreted. Note that every
    result of the column search must pass the delimiting condition, but not
    necessarily atomic or permissive condition.

    """
    if fs is None:
        fs = set().union(*[set().union(*rxn[:2]) for rxn in fcrn])
    if intrp is None:
        intrp = dict()
    log.info(f'Searching column at {depth=}'.center(40-(2*depth), '*'))
    log.debug(f'Partial interpretation {inter_list(intrp)=}')

    sicrn = subst(icrn, intrp)
    T = updateT(fcrn, sicrn, fs)
    if not checkT(T):
        #log.debug(f'Delimiting condition not satisfied with {inter_list(intrp)=}')
        raise SpeciesAssignmentError('Delimiting condition not satisfied.')

    if unknown is None:
        unknown = [i for i in range(len(fcrn))]

    # Start with the unknown column with the minimal number of True's.
    # I.e. the implementation rxn that can be assigned to the least formal reactions.
    unknown = sorted(unknown, key = lambda x: sum(T[j][x] for j in range(len(sicrn))))

    if len(unknown) == 0: 
        # Note, we cannot test for the atomic condition at this point, because
        # some trivial reactions may be needed to satisfy the atomic condition.
        yield intrp
        return

    later = dict() # less promising partial interpretations ...
    def search(later, n = 0, m = None):
        # use this to define a search hierarchy when sending results into the
        # next round. Currently, n and m are start and stop in a range of 
        # how many species interpet to a single species.
        while later:
            if m and n > m:
                break
            if n in later:
                for (inext, unext, depth) in later[n]:
                    try:
                        for isuccess in search_column(fcrn, icrn, fs, 
                                                      inext, unext, depth + 1):
                            yield isuccess
                    except SpeciesAssignmentError:
                        continue
                del later[n]
            n += 1

    for c in unknown:
        frxn = fcrn[c]
        unext = [u for u in unknown if u != c]
        for k, irxn in enumerate(sicrn):
            if not T[k][c]:
                continue
            log.debug(f'Interpret: {pretty_rxn(rl(irxn))} => {pretty_rxn(rl(frxn))}')
            # left 
            ul = irxn[0] - frxn[0]
            sl = frxn[0] - irxn[0]
            # right
            ur = irxn[1] - frxn[1]
            sr = frxn[1] - irxn[1]

            if len(ul) or len(ur):
                [kl, vl] = zip(*ul.items()) if len(ul) else [[], []]
                [kr, vr] = zip(*ur.items()) if len(ur) else [[], []]
                tmpl = enumL(len(ul), list(sl.elements()), vl)
                tmpr = enumL(len(ur), list(sr.elements()), vr)
                
                for i, j in product(tmpl, tmpr):
                    intrpleft = dict(zip(kl, i))
                    intrpright = dict(zip(kr, j))
                    for key in set(intrpleft) & set(intrpright):
                        if any([sorted(intrpleft[key][fsp]) != \
                                sorted(intrpright[key][fsp]) for fsp in fs]):
                            # Incompatible dictionaries!
                            break
                    else:
                        inext = intrp.copy()
                        inext.update({k: Counter(v) for k, v in intrpleft.items()})
                        inext.update({k: Counter(v) for k, v in intrpright.items()})
                        level = max(len(v) for k, v in chain(intrpleft.items(), 
                                                             intrpright.items()))
                        later[level] = (later.get(level, [])) + [(inext, unext, depth)]
            else: # Just in case the full interpretation is provided.
                later[0] = (later.get(0, [])) + [(intrp.copy(), unext, depth)]
            for isuccess in search(later, n = 0, m = 1):
                yield isuccess
        for isuccess in search(later, n = 1, m = 2):
            yield isuccess
    for isuccess in search(later, n = 2):
        yield isuccess
    assert len(later) == 0
    return

def crn_bisimulations(fcrn, icrn, 
                      interpretation = None,
                      formals = None, 
                      searchmode = 'time-efficient',
                      permissive = 'graphsearch',
                      permissive_depth = None):
    """ Iterate over all crn bisimulations.

    for e, bisim in enumerate(crn_bisimulations(fcrn, icrn), 1):
        print(f'Bisimulation {e} = {bisim}')

    Arguments:
        fcrn (list): A formal CRN.
        icrn (list): An implementation CRN.
        interpretation (dict, optional): A (partial) interpretation of 
            the format interpretation[is] = list[fs,..].
            Defaults to None. 
        formals (set, optional): The set of formal species. Defaults to all
            species in the formal CRN.
        permissive (string, optional): Select the method to check the
            permissive condition:
             - 'graphsearch': construct a reachability graph for each formal
                reaction.  Uses poly(n^k) space and time, where n is size of
                CRN and k is number of reactants in formal reaction.
             - 'loopsearch': search for productive loops with space-efficient
                algorithm.  Uses poly(n,k) space and poly(2^n,2^k) time.
             - 'reactionsearch': depth-first search for a path to implement each
                formal reaction.  Space and time bounds not known, but probably worse.
            Defaults to graphsearch.
        permissive_depth (int, optional): A bound on a quantity which is
            approximately the length of a path to search for, but depends on
            which algorithm is used. Defaults to None.

    Outputs:
        Yields all correct CRN bisimulations, or None if no CRN bisimulation exists.
    """
    if searchmode not in ('time-efficient', 'space-efficient'):
        raise CRNBisimulationError(f'Unsupported CRN-bisimulation search mode: {searchmode}')
    if interpretation is None:
        interpretation = dict()
    if formals is None:
        formals = set().union(*[set().union(*rxn[:2]) for rxn in fcrn])

    log.debug('Testing:')
    log.debug('Original formal CRN:')
    [log.debug('  {}'.format(pretty_rxn(r))) for r in fcrn]
    log.debug('Original implementation CRN:')
    [log.debug('  {}'.format(pretty_rxn(r))) for r in icrn]
    log.debug(f'Formal species: {formals}')

    if icrn == []: # Empty implementation CRN!
        yield {} if fcrn == [] else [{}, 0, [[], []]]

    if permissive not in ['graphsearch', 'loopsearch', 'reactionsearch']:
        raise CRNBisimulationError('Uknown option: {}'.format(
            'the permissive test should be {}'.format(
            '"graphsearch", "loopsearch", or "reactionsearch".')))
    new = []
    for [r, p] in icrn:
        nr = [deformalize(k, formals) for k in r]
        np = [deformalize(k, formals) for k in p]
        new.append([nr, np])
    icrn = new
    log.debug('Internal implementation CRN:')
    [log.debug('  {}'.format(pretty_rxn(r))) for r in icrn]
    inter = {deformalize(k, formals): v for k, v in interpretation.items()
                                                } if interpretation else {}
    log.debug('Internal interpretation:')
    [log.debug('  {}: {}'.format(k,v)) for (k, v) in inter.items()]

    if permissive_depth:
        permissive = [permissive, permissive_depth]
    log.debug(f'Permissive argument: {permissive}')

    # Last point to move to internal data structures ...
    fcrn = [[Counter(part) for part in rxn] for rxn in fcrn]
    icrn = [[Counter(part) for part in rxn] for rxn in icrn]
    intrp = inter_counter(inter)

    log.info(f'Searching for bisimulation.')
    for parti in search_column(fcrn, icrn, formals, intrp):
        try:
            for bisim in search_row(fcrn, icrn, formals, parti, 
                                    mode = searchmode,
                                    permcheck = permissive):
                yield inter_list(formalize(bisim))
        except SpeciesAssignmentError:
            continue
    return

def modular_crn_bisimulation_test(fcrns, icrns, formals, 
                                  interpretation = None, 
                                  searchmode = 'time-efficient',
                                  permissive = 'graphsearch',
                                  permissive_depth = None):
    """ Check if a modulular CRN bisimulation exists. 

    Note: There are a few modifications to the original source:
        - the arguments to check CRN bisimulation for each module changed:
            - only the formal species present in the module are passed on
        - the modularity condition input changed
            - isc (former ispCommon) are now all implementation species that
              are both in the current module and in at least one other module.
            - fsc are all formal species that are both in the current module
              and in at least one other module.
    Args:
        same as for crn_bisimulations()

    Output:
        [bool, dict]: Whether a modular CRN bisimulation exists, and the interpretation.
    """
    # Let's make a copy here to avoid modification of the partial interpretation.
    inter = {k : v for k, v in interpretation.items()} if interpretation else {}

    # Identify common implementation species.
    ispc = dict() # Store for every implementation species in which module it appears: 
    for e, module in enumerate(icrns, 1):
        mspecies = set().union(*[set().union(*rxn[:2]) for rxn in module])
        for isp in mspecies:
            ispc[isp] = ispc.get(isp, []) + [e]
    log.debug(f'ispc = {ispc}')

    # Identify common formal species.
    fspc = dict() # Store for every formal species in which module it appears: 
    for e, module in enumerate(fcrns, 1):
        mspecies = set().union(*[set().union(*rxn[:2]) for rxn in module])
        for fsp in mspecies:
            fspc[fsp] = fspc.get(fsp, []) + [e]
    log.debug(f'fspc = {fspc}')

    for e, (fcrn, icrn) in enumerate(zip(fcrns, icrns), 1):
        # Prepare inputs for crn bisimulation of this module
        mfs = {k for k in formals if e in fspc[k]}
        minter = {k: v for k, v in inter.items() if e in ispc[k]}
        found = False
        for bisim in crn_bisimulations(fcrn, icrn, 
                                       interpretation = minter, 
                                       formals = mfs, 
                                       searchmode = searchmode,
                                       permissive = permissive, 
                                       permissive_depth = permissive_depth):
            # Get all formal and implementation species that are in
            # common with at least one other module.
            fsc = {f for f, m in fspc.items() if e in m and len(m) > 1}
            isc = {i for i, m in ispc.items() if e in m and len(m) > 1}
            if passes_modularity_condition(bisim, icrn, isc, fsc):
                found = True
                inter.update(bisim)
                break # TODO: maybe re-insert the iterate part here ...
            log.debug(f'Skipping non-modular bisimulation: {bisim}')
        if not found:
            return False, None
    return True, inter

def crn_bisimulation_test(fcrn, icrn, formals, 
                          interpretation = None,
                          searchmode = 'time-efficient',
                          permissive = 'graphsearch',
                          permissive_depth = None):
    """ Backward compatible CRN bisimulation interface.

    Args:
        same as for crn_bisimulations()

    Output:
        [bool, dict]: Whether a modular CRN bisimulation exists, and the interpretation.
    """
    iterator = crn_bisimulations(fcrn, icrn, 
                                 interpretation = interpretation,
                                 formals = formals, 
                                 searchmode = searchmode,
                                 permissive = permissive,
                                 permissive_depth = permissive_depth)

    try:
        bisim = next(iterator)
        return True, bisim
    except StopIteration:
        return False, None

