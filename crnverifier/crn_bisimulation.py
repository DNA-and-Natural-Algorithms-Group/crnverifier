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
import random
from collections import Counter
from itertools import product, combinations, chain

from .utils import pretty_rxn
from .utils import interpret as interpretL

class SpeciesAssignmentError(Exception):
    pass

class SearchDepthExceeded(Exception):
    pass

class CRNBisimulationError(Exception):
    pass

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

# Utils for list based multiset operations -> it's faster than counters ...
def subsetsL(x, reverse = True):
    """ Generate all (uniqe) multi-subsets of a list x. """
    if reverse:
        return chain(*[combinations(x, l) for l in range(len(x), -1, -1)])
    else:
        return chain(*[combinations(x, l) for l in range(0, len(x) + 1)])

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

def subst(crn, intrp):
    """ Substitute implementation species with formal species in CRN.  """
    return [[interpretL(part, intrp) for part in rxn] for rxn in crn]

def same_reaction(irxn, frxn, fs):
    """ Check if irxn *could* implement frxn.

    This assumes that irxn is already interpreted using the formal species fs.
    SB: Note, that I introduced a new expression expression here, see below.
    """
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

def makeT(fcrn, icrn, fs):
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
        rr = [same_reaction(irxn, frxn, fs) for frxn in fcrn]
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

def minimal_implementation_states(fstate, inter):
    """ Generate set of minimal implementation states for a formal state.

    The returned implementation states may interpret to supersets of fstate.
    """
    def gen_supersets():
        if len(fstate) == 0:
            yield []
        else:
            fs = set(fstate).pop()
            fs1 = fstate[0]
            for k in inter:
                if fs1 in inter[k]:
                    for out in minimal_implementation_states(subtractL(fstate, 
                                                                       inter[k], 
                                                                       strict = False), inter):
                        yield [k] + out

    def filter_minimals(lst):
        seen = set()
        for x in lst:
            z = tuple(sorted(x))
            if z in seen:
                continue
            if len(x) == 0:
                yield x
                seen.add(z)
                continue
            for y in combinations(x, len(x)-1):
                if is_contained(fstate, interpretL(y, inter)):
                    break
            else:
                yield x
            seen.add(z)
        return
    return filter_minimals(gen_supersets())

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

    trivial_irxns = [] # trivial reactions
    for rxn in icrn:
        iR = interpretL(rxn[0], bisim) 
        iP = interpretL(rxn[1], bisim) 
        if sorted(iR) == sorted(iP):
            trivial_irxns.append(rxn)

    done = {s for s in bisim if Y(s) or Z(s)}
    todo = {s: [set(), set()] for s in bisim if s not in done} 
    # todo[s] = [reach, produce]

    reset = len(todo) > 0
    while reset:
        reset = False
        for r in list(todo.keys()): # start at any of the non-common intermediates
            [r_reach, r_produce] = todo[r]
            for R, P in trivial_irxns:
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

def check_permissive_bruteforce(fcrn, icrn, fs, inter, maxdepth = 100):
    """ Check the permissive condition via a brute force path search.
    
    Args:
        fcrn: The formal CRN.
        icrn: The implementation CRN.
        fs: The formal species.
        inter: The interpretation.
        maxdepth: A maximum brute force recursion depth. Defaults to 100.

    Raises:
        SearchDepthExceeded: If the hardcoded maxdepth is reached.

    """
    sicrn = subst(icrn, inter)
    T = makeT(fcrn, sicrn, fs)
    assert checkT(T) # Assuming this has been checked before calling permissive.
    assert all(fs | set(v) == fs for k, v in inter.items())

    log.debug(f'The implementation CRN:\n' + '\n'.join(
        [f'{e} {t=} {pretty_rxn(sicrn[e])} ({pretty_rxn(icrn[e])})' \
                for e, t in enumerate(T)]))
    log.debug(f'{inter=}')

    def bruteforce(si, seen, trivials, depth = 0):
        if maxdepth and depth > maxdepth:
            raise SearchDepthExceeded('WARNING: brute-force permissive test could not find a ')
        st = tuple(sorted(si))
        if st in seen:
            return False
        seen.add(st)
        for irxn in formal_irxns:
            if len(subtractL(irxn[0], si, False)) == 0:
                return True
        for irxn in trivials:
            if len(subtractL(irxn[0], si, False)) == 0:
                nS = subtractL(si, irxn[0]) + irxn[1]
                # Shuffle the crn to reduce the probability of getting stuck in
                # an infinite state space assuming that there is a way out ...
                nexttrivials = trivials[:]
                random.shuffle(nexttrivials)
                if bruteforce(nS, seen, nexttrivials, depth + 1):
                    return True
        return False

    # Find implementation reactions that interpret to trivial reactions.
    trivial_irxns = [irxn for e, irxn in enumerate(icrn) if T[e][-1]]
    for i, frxn in enumerate(fcrn):
        # Find implementation reactions that interpret to this formal reaction.
        formal_irxns = [irxn for e, irxn in enumerate(icrn) if T[e][i]]
        for S0 in minimal_implementation_states(frxn[0], inter):
            assert is_contained(frxn[0], interpretL(S0, inter))
            found = bruteforce(S0, set(), trivial_irxns)
            if not found: 
                return False, (frxn, S0)
    return True, maxdepth

def check_permissive_loopsearch(fcrn, icrn, fs, inter):
    """ The 'loopsearch' algorithm to check the permissive condition. 

    TODO: Acutally this is a very different algorithm then described in the
    paper, so before using this in production, one has to go back and see that
    it is correct and what the time/space complexity of this algorithm is.
    """
    log.warning('This "loopsearch" algorithm is not implemented as described in JDW (2019).')
    sicrn = subst(icrn, inter)
    T = makeT(fcrn, sicrn, fs)
    assert checkT(T) # Assuming this has been checked before calling permissive.
    assert all(fs | set(v) == fs for k, v in inter.items())

    log.debug(f'The implementation CRN:\n' + '\n'.join(
        [f'{e} {t=} {pretty_rxn(sicrn[e])} ({pretty_rxn(icrn[e])})' \
                for e, t in enumerate(T)]))
    log.debug(f'{inter=}')

    # Potential null species.
    nulls = [k for k, v in inter.items() if len(v) == 0]
    # Find implementation reactions that interpret to trivial reactions.
    trivial_irxns = [irxn for e, irxn in enumerate(icrn) if T[e][-1]]

    def reach_with_inf(start, goal = None, pickup = None, ignore = None, k = 0):
        """
        Search for a path from start to goal (states of non-null species) which
        produces at least pickup null species assuming it already has infinite
        copies of everything in ignore of length at most 2^k if goal is None,
        the goal is any reaction in formal_irxns.
        
        Args:
            start: an implementation state corresponding to a formal state.
            goal: a (different) implementation state corresponding to the same
                formal state, or None. If goal is None, then the goal is to find
                any implementation reaction corresponding to the formal reaction.
            pickup: xxx Something related to null species.
            ignore: Assuming there exist infinite copies of those species.
            k: Maximum search depth.
        """
        log.debug(f'CHECK: {start=} {goal=} {ignore=} {pickup=} {k=}')
        if ignore is None:
            ignore = set()
        if pickup is None:
            pickup = []

        #log.debug(f'{start=} {goal=} {ignore=}, {pickup=}, {k=}')
        if goal is None: 
            # Check if a implementation reaction is possible.
            for irxn in formal_irxns:
                # if irxn can happen: 
                if ignore >= set(subtractL(irxn[0], start, False)):
                    log.debug(f'True formal {irxn=}')
                    return True
        elif k == 0:
            # See if start can reach a pickup goal by trivial reactions.
            for irxn in trivial_irxns:
                # if irxn can happen: 
                if ignore >= set(subtractL(irxn[0], start, False)):
                    # and if goal is contained in the results of the trivial rxn.
                    if is_contained(goal+pickup, subtractL(start, irxn[0], False) + irxn[1]):
                        return True
        else:
            if reach_with_inf(start, goal, pickup, ignore, k-1):
                # If we can reach goal with trivial reactions, do it!
                log.debug(f'Shortcut {start=} {goal=} {pickup=} {ignore=} {k=}')
                return True

            # Divide the search into two by choosing a new state mid: 
            #   - find path from start to mid.
            #   - find path from mid to goal.
            for pickpart in map(list, subsetsL(pickup)):
              for mid in M:
                if mid == start or mid == goal: 
                    continue
                log.debug(f'Reach w inf: {start=} {mid=} {goal=} {pickpart=} {k=}')
                if reach_with_inf(start, mid, pickpart, ignore, k-1):
                    if reach_with_inf(mid, goal, subtractL(pickup, pickpart), ignore, k-1):
                        return True
        log.debug(f'FALSE: {start=} {goal=} {ignore=} {pickup=} {k=}')
        return False

    for i, frxn in enumerate(fcrn):
        # Find implementation reactions that interpret to this formal reaction.
        formal_irxns = [irxn for e, irxn in enumerate(icrn) if T[e][i]]
        # From the paper M(r) is the set of minimal implementation states for a formal reaction r.
        M = list(map(sorted, minimal_implementation_states(frxn[0], inter)))

        k, roundequiv = 0, 1
        for nequiv in range(1, len(M)+1):
            if nequiv > roundequiv:
                k += 1
                roundequiv *= 2
        log.debug(f'LS: {k=} to find paths between {M} states: {frxn=}')

        for S0 in M:
            # The formal reaction must be able to happen.
            assert is_contained(frxn[0], interpretL(S0, inter))

            # Is the formal reaction possible w/o null species? 
            if reach_with_inf(S0, None):
                # This state can exit by itself.
                continue

            found = False
            reset = True
            ignore = set()
            # Can we reach a state from where the formal reaction is possible?
            while reset:
                log.debug(f'{reset=}, {found=}, {ignore=}')
                reset = False

                for pickup in map(list, nulls):
                    if set(pickup) & ignore:
                        continue
                    assert len(pickup) != 0
                    log.debug(f'{S0=}, {pickup=}, {ignore=}, {k=} ')

                    if reach_with_inf(S0, S0, pickup, ignore, k):
                        log.debug(f'Pickup! {S0=} {pickup=} {ignore=}.')
                        ignore |= set(pickup)
                        reset = True
                        if reach_with_inf(S0, None, ignore = ignore):
                            found = True
                            break

                for Si in M:
                    if Si == S0: 
                        continue
                    log.debug(f'Testing: {S0=} {Si=} {ignore=} {k=}.')
                    if reach_with_inf(S0, Si, ignore = ignore, k = k):
                        log.debug(f'True: {S0=} {Si=} {ignore=} {k=}.')
                        if reach_with_inf(Si, None, ignore = ignore):
                            log.debug(f'Found exit: {Si=}.')
                            found = True
                            break
                if found is True:
                    break

            if not found:
                return False, S0
    return True, None

def check_permissive_graphsearch(fcrn, icrn, fs, inter):
    """ The 'graphsearch' algorithm to check the permissive condition. 

    Checks whether every implementation of a formal state that can react in the
    formal CRN by reaction r can also react in the implementation CRN with a 
    reaction that interprets to r.
    """

    sicrn = subst(icrn, inter)
    T = makeT(fcrn, sicrn, fs)
    assert checkT(T) # assuming this has been checked before calling permissive.
    assert all(fs | set(v) == fs for k, v in inter.items())

    log.debug(f'The implementation CRN:\n' + '\n'.join(
        [f'{e} {t=} {pretty_rxn(sicrn[e])} ({pretty_rxn(icrn[e])})' \
                for e, t in enumerate(T)]))
    log.debug(f'{inter=}')

    # Species that interpret to nothing.
    nullsp = set([k for k, v in inter.items() if len(v) == 0])

    # Find implementation reactions that interpret to trivial reactions.
    trivial_irxns = [irxn for e, irxn in enumerate(icrn) if T[e][-1]]

    max_depth = 0
    for i, frxn in enumerate(fcrn):
        # Find implementation reactions that interpret to this formal reaction.
        formal_irxns = [irxn for e, irxn in enumerate(icrn) if T[e][i]]
        # All implementation states which must permit the formal reaction.
        M = list(map(tuple, map(sorted, minimal_implementation_states(frxn[0], inter))))
        done = set() # set of species known to implement the current formal rxn.
        todo = {istate: [set(), set()] for istate in M}
        changed, depth = True, 0
        while changed:
            changed, depth = False, depth + 1
            for istate in M:
                if istate in done:
                    continue
                for irxn in formal_irxns:
                    # If the implementation state can implement the formal reaction.
                    if set(subtractL(irxn[0], istate, False)) <= todo[istate][0]: 
                        done.add(istate)
                        del todo[istate]
                        changed = True
                        break
                if istate in done:
                    continue
                for irxn in trivial_irxns:
                    # If the implementation state allows the implementation reaction.
                    if set(subtractL(irxn[0], istate, False)) <= todo[istate][0]: 
                        # apply the reaction to the current state.
                        nstate = subtractL(list(istate), irxn[0], False) + irxn[1]
                        for istate2 in M:
                            # if one of the other implementation states is reachable
                            if is_contained(istate2, nstate):
                                if istate2 in done:
                                    done.add(istate)
                                    del todo[istate]
                                    changed = True
                                    break
                                if istate in todo[istate2][1]: # I assume its a loop?
                                    s = todo[istate2][0] | (set(nstate) & nullsp)
                                    if not s <= todo[istate][0]:
                                        todo[istate][0] |= s
                                        changed = True
                                if istate2 not in todo[istate][1]:
                                    todo[istate][1].add(istate2)
                                    changed = True
                                if not (todo[istate2][1] <= todo[istate][1]):
                                    todo[istate][1] |= todo[istate2][1]
                                    changed = True
                    if istate in done:
                        break
        if todo: 
            return False, (frxn, M)
        if max_depth < depth: 
            max_depth = depth
    return True, max_depth
        
def passes_permissive_condition(fcrn, icrn, fs, inter, permcheck = 'graphsearch'):
    """
    The difficulty of checking the permissive condition scales with the number
    of minimal states for any given formal reaction r = (R, P), which typically
    scales like (and scales no worse than) n**k where n = |S'| (the species of
    the implementation CRN) and k = |R|. 
    
    (There may be order n**k minimal implementation states for r and the
    trajectories by which any one implements r may have to pass through most or
    all of them.  When k is unbounded, checking an interpretation is
    PSPACE-complete.)

    We provide three algorithms: 
        * The loopsearch algorithm, which runs in poly(n, k) space and poly(n**(kn)) time.
        * The graphsearch algorithm which takes poly(n**k) space and time. 

    In practice, many CRNs have large numbers of species but small numbers of
    reactants per (formal) reaction. For those CRNs, we recommend the
    graphsearch algorithm which is faster than the loopsearch algorithm when k
    is small but taking much more space when k is large.
    """
    log.debug(f'Checking permissive condition using {permcheck=}.')
    if permcheck == 'graphsearch':
        passes, info = check_permissive_graphsearch(fcrn, icrn, fs, inter)
    elif permcheck == 'loopsearch':
        passes, info = check_permissive_loopsearch(fcrn, icrn, fs, inter)
    elif permcheck == 'bruteforce':
        passes, info = check_permissive_bruteforce(fcrn, icrn, fs, inter)
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
    T = makeT(fcrn, sicrn, fs)
    return checkT(T)

def passes_atomic_condition(inter, fs):
    """ Tests an interpretation for the atomic condition.

    For every formal species there exists an implementation species which
    interprets to it.  
    """
    return all(any([f] == v for v in inter.values()) for f in fs)

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

    log.debug(f'Searching row at {depth=}'.center(80-(2*depth), '~'))
    log.debug(f'Partial interpretation {intrp=}')
    log.debug(f'Forbidden moves {moves=}')

    sicrn = subst(icrn, intrp)
    T = makeT(fcrn, sicrn, fs)
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
        [f'{e} {t=} {pretty_rxn(sicrn[e])} ({pretty_rxn(icrn[e])})' \
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
            log.debug(f'Interpret: {pretty_rxn(irxn)} => trivial')
            fl = [x for x in subtractL(irxn[0], irxn[1], False) if x in fs]
            fr = [x for x in subtractL(irxn[1], irxn[0], False) if x in fs]
            ul = Counter([x for x in subtractL(irxn[0], irxn[1], False) if x not in fs])
            ur = Counter([x for x in subtractL(irxn[1], irxn[0], False) if x not in fs])
            if len(ul) or len(ur):
                [kl, vl] = zip(*ul.items()) if len(ul) else [[], []]
                tmpl = enumL(len(ul), fr, vl)
                [kr, vr] = zip(*ur.items()) if len(ur) else [[], []]
                tmpr = enumL(len(ur), fl, vr)

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
                        inext.update({k: list(v) for k, v in imove})
                        level = max(len(v) for k, v in imove)
                        if level in later:
                            later[level].append((inext, depth, imove))
                        else:
                            later[level] = [(inext, depth, imove)]

        for c, frxn in enumerate(fcrn): # Assign a formal reaction.
            if not T[k][c]:
                continue
            log.debug(f'Interpret: {pretty_rxn(irxn)} => {pretty_rxn(frxn)}')
            # left 
            ul = Counter(subtractL(irxn[0], frxn[0], False))
            sl = subtractL(frxn[0], irxn[0], False)
            # right
            ur = Counter(subtractL(irxn[1], frxn[1], False))
            sr = subtractL(frxn[1], irxn[1], False)

            if len(ul) or len(ur):
                [kl, vl] = zip(*ul.items()) if len(ul) else [[], []]
                [kr, vr] = zip(*ur.items()) if len(ur) else [[], []]
                tmpl = enumL(len(ul), sl, vl)
                tmpr = enumL(len(ur), sr, vr)

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
                        inext.update({k: list(v) for k, v in imove})
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

    This "column search" finds all combinations where all formal reactions have
    one matching implementation reaction. For example, if there is one formal
    reaction and three compatible implementation reactions, then this function
    returns three interpretations that have exactly one reaction interpreted.
    Note that every result of the column search must pass the delimiting
    condition, but not necessarily atomic or permissive condition.
    """
    if fs is None:
        fs = set().union(*[set().union(*rxn[:2]) for rxn in fcrn])
    if intrp is None:
        intrp = dict()
    log.debug(f'Searching column at {depth=}'.center(40-(2*depth), '*'))
    log.debug(f'Partial interpretation {intrp=}')

    sicrn = subst(icrn, intrp)
    T = makeT(fcrn, sicrn, fs)
    if not checkT(T):
        #log.debug(f'Delimiting condition not satisfied with {intrp=}')
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
            log.debug(f'Interpret: {pretty_rxn(irxn)} => {pretty_rxn(frxn)}')
            # left 
            ul = Counter(subtractL(irxn[0], frxn[0], False))
            sl = subtractL(frxn[0], irxn[0], False)
            # right
            ur = Counter(subtractL(irxn[1], frxn[1], False))
            sr = subtractL(frxn[1], irxn[1], False)

            if len(ul) or len(ur):
                [kl, vl] = zip(*ul.items()) if len(ul) else [[], []]
                [kr, vr] = zip(*ur.items()) if len(ur) else [[], []]
                tmpl = enumL(len(ul), sl, vl)
                tmpr = enumL(len(ur), sr, vr)
                
                for i, j in product(tmpl, tmpr):
                    intrpleft = {k: tuple(sorted(v)) for k, v in zip(kl, i)}
                    intrpright = {k: tuple(sorted(v)) for k, v in zip(kr, j)}
                    for key in set(intrpleft) & set(intrpright):
                        if any([sorted(intrpleft[key][fsp]) != \
                                sorted(intrpright[key][fsp]) for fsp in fs]):
                            # Incompatible dictionaries!
                            break
                    else:
                        inext = intrp.copy()
                        inext.update({k: list(v) for k, v in intrpleft.items()})
                        inext.update({k: list(v) for k, v in intrpright.items()})
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
                      permissive = 'graphsearch'):
    """ Iterate over all crn bisimulations.

    for e, bisim in enumerate(crn_bisimulations(fcrn, icrn), 1):
        print(f'Bisimulation {e} = {bisim}')

    Args:
        fcrn (list): A formal CRN.
        icrn (list): An implementation CRN.
        interpretation (dict, optional): A (partial) interpretation with the
            format interpretation[is] = list[fs,..]. Defaults to None. 
        formals (set, optional): The set of formal species. 
            Defaults to None: all species in the formal CRN.
        permissive (string, optional): Select the method to check the
            permissive condition:
             - 'graphsearch': construct a reachability graph for each formal
               reaction.  Uses poly(n^k) space and time, where n is size of
               implementation CRN and k is the maximum number of reactants in
               any formal reaction.
             - 'loopsearch': search for productive loops with space-efficient
                algorithm. Uses poly(n, k) space and poly(2^n, 2^k) time.
             - 'bruteforce': a depth-first search for a path to implement each
               formal reaction. Space and time bounds are not known, but
               probably worse than 'graphsearch' and 'loopsearch'.
            Defaults to None: 'graphsearch' or 'loopsearch' dependent on the
                values of n and k.
    Yields:
        All correct interpretations (CRN bisimulations).
    """
    if searchmode not in ('time-efficient', 'space-efficient'):
        raise CRNBisimulationError(f'Unsupported CRN-bisimulation search mode: {searchmode}')
    if interpretation is None:
        interpretation = dict()
    if formals is None:
        formals = set().union(*[set().union(*rxn[:2]) for rxn in fcrn])

    log.debug(f'Finding CRN bisimulations using {permissive} algorithm:')
    log.debug('Original formal CRN:')
    [log.debug('  {}'.format(pretty_rxn(r))) for r in fcrn]
    log.debug('Original implementation CRN:')
    [log.debug('  {}'.format(pretty_rxn(r))) for r in icrn]
    log.debug(f'Formal species: {formals}')

    if icrn == []: # Empty implementation CRN!
        yield {} if fcrn == [] else [{}, 0, [[], []]]

    if permissive not in ['graphsearch', 'loopsearch', 'bruteforce']:
        raise CRNBisimulationError('Unknown option: {}'.format(
            'the permissive test should be {}'.format(
            '"graphsearch", "loopsearch", or "bruteforce".')))
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


    log.info(f'Searching for bisimulation.')
    for parti in search_column(fcrn, icrn, formals, inter):
        try:
            for bisim in search_row(fcrn, icrn, formals, parti, 
                                    mode = searchmode,
                                    permcheck = permissive):
                yield formalize(bisim)
        except SpeciesAssignmentError:
            continue
    return

def modular_crn_bisimulation_test(fcrns, icrns, formals, 
                                  interpretation = None, 
                                  permissive = 'graphsearch'):
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
                                       permissive = permissive):
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
                          permissive = 'graphsearch'):
    """ Backward compatible CRN bisimulation interface.

    Args:
        same as for crn_bisimulations()

    Output:
        [bool, dict]: Whether a modular CRN bisimulation exists, and the interpretation.
    """
    iterator = crn_bisimulations(fcrn, icrn, 
                                 interpretation = interpretation,
                                 formals = formals, 
                                 permissive = permissive)

    try:
        bisim = next(iterator)
        return True, bisim
    except StopIteration:
        return False, None

