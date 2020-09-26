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
from itertools import product, permutations, combinations, chain

from .utils import interpret as interpretL
from collections import Counter # i.e. multiset (states of a CRN, interpretations, ...)

from .deprecated import subsets, enum, searchr

class SpeciesAssignmentError(Exception):
    pass

class CRNBisimulationError(Exception):
    pass

# I/O utils
def pretty_crn(crn):
    for rxn in crn:
        yield '{} -> {}'.format(' + '.join(rxn[0]), ' + '.join(rxn[1]))

def pretty_rxn(rxn, flag_internal = True):
    R = list(map(lambda r: 'i{{{}}}'.format(r[1]) if isinstance(
        r, tuple) else '{}'.format(r), rxn[0]))
    P = list(map(lambda r: 'i{{{}}}'.format(r[1]) if isinstance(
        r, tuple) else '{}'.format(r), rxn[1]))
    return '{} -> {}'.format(' + '.join(R), ' + '.join(P))

def inter_counter(inter):
    # For internal representation of interpretations 
    # return interpretations of the format: inter[str] = Counter
    return {k: Counter(v) for k, v in inter.items()}

def inter_list(inter):
    # For consistency with the rest of the package, 
    # return interpretations of the format: inter[str] = list
    return {k: list(v.elements()) for k, v in inter.items()}

def deformalize(k, fs):
    # Do not confuse formal species with same-named species in the
    # implementation CRN.
    return ('impl', k) if k in fs else k

def formalize(intrp):
    intr = {}
    for sp in intrp:
        if isinstance(sp, tuple):
            assert sp[0] == 'impl'
            intr[sp[1]] = intrp[sp]
        else:
            intr[sp] = intrp[sp]
    return intr

def subsetsL(x):
    """ Generate all (uniqe) multi-subsets of a list x. """
    return chain(*[combinations(x, l) for l in range(len(x), -1, -1)])

def substractL(l1, l2):
    l2 = l2[:]
    for el in l1:
        l2.remove(el)
    return l2

def is_contained(a, b):
        # True if multiset a is contained in multiset b.
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
            for j in enumL(n-1, substractL(ss, l), weights[1:]):
                yield [wss] + j
        except SpeciesAssignmentError:
            continue
    return

def solve_contejean_devie(a):
    """ Algorithm from Contejean & Devie 1994.

    Find a non-negative and non-trivial integer solution x of the equation ax=0.
    Return [] when there is no such solution.
    """
    q = len(a[0])
    
    def multi(x, y):
        s = 0
        for i in range(len(x)):
            s = s + x[i] * y[i]
        return s
    
    def sub(x):
        s = []
        for i in range(len(a)):
            s.append(multi(a[i],x))
        return s
    
    def Min(b,t):
        if not b:
            return True
        else:
            for i in range(len(b)):
                r = True;
                for j in range(q):
                    r = r and (b[i][j] <= t[j])
                if r:
                    return False
            return True
        
    e = []
    for i in range(q):
        e.append([])
        for j in range(len(a)):
            e[i].append(a[j][i])     
    p = []
    frozen = []
    for i in range(q):
        p.append([1 if j == i else 0 for j in range(q)])
        frozen.append([i == q-1 or j < i for j in range(q)])
    zero = [0 for i in range(len(a))]
    zero1 = [0 for i in range(q)]
    b = []
    while p:
        t = p.pop()
        if sub(t) == zero:
            if t[q-1] == 1:
                return t    # just get the first solution, not all solutions (unlike C&D 1994).
            b.append(list(t))
            frozen.pop()
        else:
            f = frozen.pop()
            for i in range(q):
                if not f[i] and (multi(sub(t), e[i]) < 0):
                    tmp = list(t)
                    tmp[i] += 1
                    if Min(b, tmp):
                        if i == q-1:
                            f[i] = True
                        p.append(tmp)
                        frozen.append(list(f))
                    f[i] = True
    return []

def solve(a):
    # wrapper method to solve a system of equations with method of choice
    return solve_contejean_devie(a)

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

def subst(crn, intrp):
    # Substitute implementation species for formal species in CRN according to
    # interpretation.
    return [[interpret(j, intrp) for j in rxn] for rxn in crn]

def check_delimiting(fcrn, icrn, fs, inter):
    # TODO: untested and not used, just to keep the idea ...
    sicrn = subst(icrn, inter)
    T = updateT(fcrn, sicrn, fs)
    return checkT(T)

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

def same_reaction(irxn, frxn, fs):
    """ Check if irxn *could* implement frxn.

    This assumes that irxn is already interpreted using the formal species fs.
    SB: Note, that I introduced a new expression expression here, see below.
    """
    fr = set(frxn[0] - irxn[0]) # exess formal reactants 
    ir = set(irxn[0] - frxn[0]) # exess implementation reactants 
    fp = set(frxn[1] - irxn[1]) # exess formal products
    ip = set(irxn[1] - frxn[1]) # exess implementation products
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

def updateT(fcrn, icrn, fs):
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

        ir = set(irxn[0] - irxn[1]) # "non-trivial" reactants
        ip = set(irxn[1] - irxn[0]) # "non-trivial" products
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

def perm(fcrn, icrn, fs, intrp, permcheck, state):
    """ Check the permissive condition.
    
    Uses the original formal CRN  and original implementation CRN
    (without substitutions).

    Args:
        fcrn: The original formal CRN
        icrn: The original implementation CRN
        fs: The formal species.

    """
    if len(permcheck) == 2:
        permissive_depth = permcheck[1]
        permcheck = permcheck[0]
    else:
        permissive_depth = None

    intr, max_depth, permissive_failure = state
    #log.info(f'Checking permissive condition using {permcheck=} {permissive_depth=}.')
    tr = []
    fr = []
    hasht = set([])
    crn2 = subst(icrn, intrp)
    T = updateT(fcrn, crn2, fs)
    if not checkT(T):
        return [False, state]

    nulls = [k for k in intrp if not len(list(intrp[k]))] # null species

    def cnstr(s):
        # given formal state s, generate minimal set of impl. states which
        #  interpret to a superset of s
        # = concatenation of any x such that s[0] in m(x) with
        #  all minimal implementations of s - m(x)
        if list(s.keys()) == []:
            yield Counter()
        else:
            s1 = list(s.keys())[0]
            for k in intrp:
                if s1 in intrp[k]:
                    if (s - intrp[k])[s1] >= s[s1]:
                        assert False
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
        for i in tr:
            if list((i[0] - s).keys()) == []:
                t = (s - i[0]) + i[1]
                out = reactionsearch(t, d+1)
                if out:
                    return True
                elif out is None:
                    ret = None
        return ret

    def midsearch(start, goal, pickup, ignore, formal, k):
        # search for a path from start to goal (states of non-null species)
        #  which produces at least pickup null species
        #  assuming it already has infinite copies of everything in ignore
        #  of length at most 2^k
        # if goal is None, the goal is any reaction in fr
        if goal is not None:
            if ignore.issuperset(goal - start):
                return True
            if not interleq(goal, start, intrp):
                return False

        if k == 0:
            if goal is None:
                for rx in fr:
                    if ignore.issuperset(rx[0] - start):
                        return True
                return False
            
            for rx in tr:
                if ignore.issuperset(rx[0] - start):
                    # every element of rx[0] - start (multiset)
                    #  is also an element of ignore (set)
                    # e.g. true on rx[0] - start = {|a,a,b|}, ignore = {a,b,c}
                    if msleq(goal, (start - rx[0]) + rx[1]):
                        return True

        else:
            if midsearch(start,goal,pickup,ignore,formal,k-1):
                return True
            for part in subsets(Counter(pickup)):
                for mid in cnstr(formal):
                    if midsearch(start,mid,set(part),ignore,formal,k-1) \
                       and ((not interleq(start, mid, intrp)) or
                            midsearch(mid,goal,pickup-set(part),ignore,
                                      formal,k-1)):
                        return True

        return False

    def loopsearch(start, formal):
        # search for a path from start to any reaction in fr
        if permissive_depth:
            rounds = math.ceil(math.log(permissive_depth, 2))
        else:
            nequiv = 0
            rounds = 0
            roundequiv = 1
            for point in cnstr(formal):
                nequiv += 1
                if nequiv > roundequiv:
                    rounds += 1
                    roundequiv *= 2

        for parti in enum(len(nulls) + 1,Counter(nulls)):
            part = list(map(set,parti))
            if any([part[i] != set() and part[i+1] == set() \
                    for i in range(len(part) - 2)]):
                continue # avoid redundancy

            pickups = [_f for _f in part[:-1] if _f is not None]
            check1 = True
            place = start
            ignore = set()
            for pickup in pickups:
                check2 = False
                for base in cnstr(formal):
                    if midsearch(place,base,set(),ignore,formal,rounds):
                        if not interleq(place,base,intrp):
                            return True
                        elif midsearch(base,base,pickup,ignore,formal,rounds):
                            check2 = True
                            place = base
                            break

                if not check2:
                    check1 = False
                    break

                ignore |= pickup

            if check1 and midsearch(place,None,set(),ignore,formal,rounds):
                return True

        return False if not permissive_depth else None

    def graphsearch(formal):
        """Check whether every implementation of state "formal" can implement the given formal reaction.

        This function stores for each state the list of states it can reach.
         -- w/o null species (except those producible by loops)
        """
        points = list([[x,set([]),[]] for x in cnstr(formal)])
        # points[i][0] is the ith implementation state
        # points[i][1] is all null species loopable from that state
        # points[i][2][j] is True if points[j][0] is known to be reachable
        #  (from state i without initial null species)
        # exception: points[i][2] is True if a formal reaction is known to
        #  be reachable from state i
        l = len(points) # let's not recalculate this every time...
        rngl = list(range(l)) # same here...
        for i in rngl:
            points[i][2] = l*[False]

        if permissive_depth:
            d = 0
        changed = True
        while changed and ((not permissive_depth) or d < permissive_depth):
            changed = False
            if permissive_depth:
                d = d + 1
            for i in rngl:
                if points[i][2] is not True:
                    for rx in fr:
                        if points[i][1].issuperset(rx[0] - points[i][0]):
                            points[i][2] = True
                            changed = True
                            break

                    if points[i][2] is True:
                        continue

                    for rx in tr:
                        if points[i][1].issuperset(rx[0] - points[i][0]):
                            left = points[i][0] - rx[0]
                            after = left + rx[1]
                            for j in rngl:
                                if msleq(points[j][0],after):
                                    if points[j][2] is True:
                                        points[i][2] = True
                                        changed = True
                                        break

                                    if points[j][2][i]:
                                        s = points[j][1].union(
                                            [x for x in left if x in nulls])
                                        if not s <= points[i][1]:
                                            points[i][1] |= s
                                            changed = True

                                    if not points[i][2][j]:
                                        points[i][2][j] = True
                                        changed = True

                                    for k in rngl:
                                        if (not points[i][2][k]) \
                                           and points[j][2][k]:
                                            points[i][2][k] = True
                                            changed = True

                        if points[i][2] is True:
                            break

        if permissive_depth and changed:
            return None
        return all([p[2] is True for p in points])
    
    n = 0
    # build tr just once
    for i in T:
        if i.pop():
            tr.append(icrn[n])
        n += 1
    for i in range(len(fcrn)):
        # build fr for this formal reaction
        fr = []
        for j in range(len(icrn)):
            if T[j][i]:
                fr.append(icrn[j])

        ist = cnstr(fcrn[i][0])

        if permcheck == "graphsearch":
            out = graphsearch(fcrn[i][0])
            if not out:
                if max_depth == -2:
                    return [False, [intr, max_depth, permissive_failure]]
                if out is False: # proven unreachable
                    max_depth = -1
                else: # arbitrarily specified max search depth reached
                    max_depth = -2
                intr = intrp.copy()
                permissive_failure[0] = fcrn[i]
                permissive_failure[1] = ["Somewhere"]
                return [False, [intr, max_depth, permissive_failure]]

            continue
       
        # At this point, "ist" contains an exhaustive and sufficient list of
        # possible initial implementation states in which the current formal
        # reaction #i must be able to fire (via a series of trivial reactions),
        # Note that we will only want to test states in "ist" that interpret to
        # a state in which #i can fire.

        tested = []  # to avoid testing a superset of some state that's already been tested
        spaceefficient = True # ... but only if we have space to store them
        for j in ist:
            tmp = interpret(j,intrp)
            if msleq(fcrn[i][0], tmp):  # only test if reactants j interpret to allow #i to fire
                t = False
                for k in tested:
                    # will be a 0-length loop if spaceefficient
                    if msleq(k, j):
                        t = True
                        break
                if t:
                    continue
                hasht = set([])
                found = ((permcheck=="loopsearch") and loopsearch(j,fcrn[i][0])) \
                        or ((permcheck=="reactionsearch") and reactionsearch(j,0))
                if not found: # didn't work, but keep track of stuff for user output
                    if max_depth == -2:
                        return [False, [intr, max_depth, permissive_failure]]
                    if found is False: # reaction proven unreachable
                        max_depth = -1
                    else: # arbitrarily specified max search depth reached
                        max_depth = -2
                    intr = intrp.copy()
                    permissive_failure[0] = fcrn[i]
                    permissive_failure[1] = j
                    return [False, [intr, max_depth, permissive_failure]]
                elif not spaceefficient:
                    tested.append(j)
    intr = intrp.copy()
    return [True, intr]

def solve_trivial(fcrn, icrn, fs, intrp, permcheck, state):
    # All unknown implementation reactions (i.e. those with some unknown
    # species) must be trivial.  Build the matrix for the the "solve" function,
    # to see whether the interpretation can be completed as required.  Also
    # note that "unknow" is not used; the unknown species are recalculated here
    # because "unknow" contains only those implementation reactions that have
    # not been solved by the row search, but other rows (i.e. implementation
    # reactions) may be implicitly solved by the partial interpretation.
    unknown = []
    ust = set([])
    sicrn = subst(icrn, intrp)
    for i in range(len(sicrn)):
        tmp = set(sicrn[i][0])-set(fs)
        tmp |= set(sicrn[i][1])-set(fs)
        ust |= tmp
        if tmp != set([]):
            unknown.append(i)  # list of implementation reactions with an unknown species
    us = list(ust)  # all species that remain unknown in current (partial) interpretation
    
    # check atomic condition
    atoms = set()
    for k in intrp:
        sps = list(intrp[k].keys())
        if len(sps) == 1 and intrp[k][sps[0]] == 1:
            atoms.add(sps[0])

    atomsleft = list(set(fs) - atoms)
    l = len(atomsleft)
    for assign in permutations(us, l): # works even if l == 0
        # each assign is a tuple of implementation species to be interpreted
        #  as exactly one formal species, matching the order of atomsleft
        # if l == 0 then it.permutations(us, l) == [()], a list with
        #  element which is the zero-length tuple,
        #  so this loop will run exactly once with no change to itmp
        itmp = dict(intrp)
        for i in range(l):
            assert assign[i] not in itmp
            itmp[assign[i]] = Counter({atomsleft[i]: 1})

        sicrntmp = subst(icrn, itmp)
        T = updateT(fcrn, sicrntmp, fs)
        if (not checkT(T)) or any([not T[i][-1] for i in unknown]):
            # either the table is bad, or a reaction we are assuming
            #  is trivial can no longer be trivial
            continue

        ustmp = [sp for sp in us if sp not in assign]
        if not ustmp:
            out = perm(fcrn, icrn, fs, itmp, permcheck, state)
            if out[0]:
                return out
            else:
                state = out[1]
                continue

        n = 0
        a = []
        for i in unknown:
            a.append([])
            for j in ustmp:
                a[n].append(sicrntmp[i][0][j]-sicrntmp[i][1][j])

            a[n].append(0)
            n += 1

        for isp in ustmp:
            itmp[isp] = Counter()

        check = True
        for fsp in fs:
            n = 0
            for j in unknown:
                a[n].pop()
                a[n].append(sicrntmp[j][0][fsp]-sicrntmp[j][1][fsp])
                n += 1

            s = solve(a)
            if s == []:
                check = False
                break
            else:
                for j in range(len(ustmp)):
                    itmp[ustmp[j]][fsp] = s[j]

        if check:
            for isp in ustmp:
                itmp[isp] = itmp[isp] + Counter()
            out = perm(fcrn, icrn, fs, itmp, permcheck, state)
            if out[0]:
                return out
            else:
                state = out[1]
                continue

    return [False, state]

def atomic_condition(inter, fs):
    """ Tests an interpretation for the atomic condition.

    For every formal species there exists an implementation species which interprets to it.
    """
    return all(any(set(f) == set(v.elements()) for v in inter.values()) for f in fs)

def search_row(fcrn, icrn, fs, 
               unknown, intrp, 
               depth, permcheck, state, nontriv = False):
    """ Row search: every implementation reaction must interpret to a formal reaction 
        (or be trivial).
    """
    #print(f'searching row at {depth=}')
    #print(f'{intrp=}')
    #print(f'{state=}')

    intr, max_depth = state[0:2]
    sicrn = subst(icrn, intrp)
    T = updateT(fcrn, sicrn, fs)
    if not checkT(T):
        log.debug(f'reset delim with {intrp=}')
        raise SpeciesAssignmentError('delimiting')

    if unknown is None:
        # It is a row indices, where unassigned species can be found ...
        unknown = [i for i, sirxn in enumerate(sicrn) if \
                len(set(sirxn[0]) - fs) or len(set(sirxn[1]) - fs)]

    # Are we done? Check the permissive condition ...
    if unknown == []:
        if not atomic_condition(intrp, fs):
            log.debug(f'reset atomic with {intrp=}')
            raise SpeciesAssignmentError('atomic')
        out = perm(fcrn, icrn, fs, intrp, permcheck, state)
        if out[0]:
            yield out[1]
            return
        else:
            log.debug(f'reset permissive with {intrp=}')
            raise SpeciesAssignmentError('permissive')
    else:
        # Start with the unknown row with the minimal number of True's (excluding trivial reactions)
        unknown = sorted(unknown, key = lambda x: T[x][:-1].count(True))
        assert True in T[unknown[0]]

    if max_depth >= 0 and depth > max_depth:
        intr = intrp.copy()
        max_depth = depth

    found = False

    # all unsearched rows could be trivial -- try it!
    if all(T[x][-1] for x in unknown) and not nontriv:  
        # This has its own permissive checker ...
        out = solve_trivial(fcrn, icrn, fs, intrp, permcheck,
                            [intr, max_depth] + state[2:])
        if out[0]:
            yield out[1]
            found = True
        nontriv = True

    for k in unknown:
        if not any(T[k][:-1]):
            # This can only be a trivial reaction ...
            continue
        unext = [u for u in unknown if u != k]

        # That is probably some form of optimization ...
        if T[k][-1] is True:  
            # This can be a trivial reaction, leave it be????
            try:
                inext = intrp.copy()
                for isuccess in search_row(fcrn, icrn, fs, 
                                           unext, inext, 
                                           depth, permcheck, 
                                           [intr, max_depth] + state[2:], 
                                           nontriv):
                    inext.update(isuccess.items())
                    yield inext
            except SpeciesAssignmentError:
                continue

        # This part cannot solve trivial reactions ... 
        for c in range(len(fcrn)):
            if not T[k][c]:
                continue
            # left 
            ul = sicrn[k][0] - fcrn[c][0]
            sl = fcrn[c][0] - sicrn[k][0]
            [kl, vl] = zip(*ul.items()) if len(ul) else [[], []]
            tmpl = enumL(len(ul), list(sl.elements()), vl)
            # right
            ur = sicrn[k][1] - fcrn[c][1]
            sr = fcrn[c][1] - sicrn[k][1]
            [kr, vr] = zip(*ur.items()) if len(ur) else [[], []]
            tmpr = enumL(len(ur), list(sr.elements()), vr)

            for i, j in product(tmpl, tmpr):
                intrpleft = dict(zip(kl, i))
                intrpright = dict(zip(kr, j))
                for key in set(intrpleft) & set(intrpright):
                    if any([sorted(intrpleft[key][fsp]) != sorted(intrpright[key][fsp]) \
                            for fsp in fs]):
                        # Incompatible dictionaries!
                        break
                else:
                    inext = intrp.copy()
                    inext.update({k: Counter(v) for k, v in intrpleft.items()})
                    inext.update({k: Counter(v) for k, v in intrpright.items()})
                    try:
                        for isuccess in search_row(fcrn, icrn, fs, 
                                                   unext, inext, 
                                                   depth + 1, permcheck,
                                                   [intr, max_depth] + state[2:]):
                            inext.update(isuccess.items())
                            yield inext
                    except SpeciesAssignmentError:
                        continue
        break # TODO: I don't understand why we want this, but it makes everyting soo fast!
    return

def search_column(fcrn, icrn, fs = None, inter = None, unknown = None, depth = 0):
    """ Find an interpretation seed that passes the permissive condition.

    TODO: need to check if it finds all of these seeds, or just some....

    Note: this does not search ever combination of how to seed an interpretatation,
    and it also does not return complete interpretations.... i wonder what this
    is doing ...
    
    search depth is not implemented.
    """
    #print(f'\nsearching {depth=}')
    if fs is None:
        fs = set().union(*[set().union(*rxn[:2]) for rxn in fcrn])
    if inter is None:
        inter = {}
    if unknown is None:
        unknown = [i for i in range(len(fcrn))]

    # A quick check if the delimiting condition is still satisfied.
    sicrn = subst(icrn, inter)
    T = updateT(fcrn, sicrn, fs)
    if not checkT(T):
        #print(f'reset with {inter=}')
        raise SpeciesAssignmentError

    # A quick check if we have found a successfull partial interpretation?
    if len(unknown) == 0: 
        # NOTE: Do not make sure atomic conditition is satisfied at this point.
        #if not atomic_condition(inter, fs):
        #    raise SpeciesAssignmentError
        #print(f'returning {inter=} at depth {depth=}')
        yield inter
        return

    # Start with the unknown column with the minimal number of True's.
    # I.e. the implementation rxn that can be assigned to the least formal reactions.
    unknown = sorted(unknown, key = lambda x: sum(T[j][x] for j in range(len(sicrn))))

    for c in unknown:
        unext = [u for u in unknown if u != c]
        for k in range(len(sicrn)):
            if not T[k][c]:
                continue
            # left 
            ul = sicrn[k][0] - fcrn[c][0]
            sl = fcrn[c][0] - sicrn[k][0]
            [kl, vl] = zip(*ul.items()) if len(ul) else [[], []]
            tmpl = enumL(len(ul), list(sl.elements()), vl)
            # right
            ur = sicrn[k][1] - fcrn[c][1]
            sr = fcrn[c][1] - sicrn[k][1]
            [kr, vr] = zip(*ur.items()) if len(ur) else [[], []]
            tmpr = enumL(len(ur), list(sr.elements()), vr)

            for i, j in product(tmpl, tmpr):
                intrpleft = dict(zip(kl, i))
                intrpright = dict(zip(kr, j))
                for key in set(intrpleft) & set(intrpright):
                    if any([sorted(intrpleft[key][fsp]) != sorted(intrpright[key][fsp]) \
                        for fsp in fs]):
                        # Incompatible dictionaries!
                        break
                else:
                    inext = inter.copy()
                    inext.update({k: Counter(v) for k, v in intrpleft.items()})
                    inext.update({k: Counter(v) for k, v in intrpright.items()})
                    try:
                        for isuccess in search_column(fcrn, icrn, fs, inext, unext, depth + 1):
                            inext.update(isuccess.items())
                            yield inext
                    except SpeciesAssignmentError:
                        continue
        break # TODO: let's break here as well, this is so weird ...
    return

def find_bisimulation(fcrn, icrn, fs, intrp, depth, permcheck, state):
    """
    fcrn, icrn, fs, unknown and intrp, permcheck remain constant, i think.

    Assumes that icrn and fcrn have diffferent species names, unless they are 
    the same specis!

    Depth will update, state will update.

    The inner icrn/fcrn/inter format layer

    """

    # push the data structure further in ...
    fcrn = [[Counter(part) for part in rxn] for rxn in fcrn]
    icrn = [[Counter(part) for part in rxn] for rxn in icrn]
    intrp = inter_counter(intrp)

    found = False # Changes to True if you found something ... but what?
    for inter in search_column(fcrn, icrn, fs, intrp):

        # Rollback to original interpretation if you reach the max_depth,
        # I assume this is a hack to do some sort of depth/breadth first search.
        intr, max_depth = state[0:2]
        if max_depth >= 0 and depth > max_depth:
            intr = intrp.copy()
            max_depth = depth
        try:
            for bisim in search_row(fcrn, icrn, fs, None, inter, depth, permcheck, 
                                     [intr, max_depth] + state[2:]):
                if found is False:
                    yield True
                    found = True
                inter.update(bisim.items())
                yield inter_list(inter)
        except SpeciesAssignmentError:
            continue

    if not found:
        yield False
        yield [intrp, 0] + state[2:]

def is_modular(bisim, icrn, common_is, common_fs):
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

def crn_bisimulations(fcrn, icrn, 
                      interpretation = None,
                      formals = None, 
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
    if interpretation is None:
        interpretation = dict()
    if formals is None:
        formals = set().union(*[set().union(*rxn[:2]) for rxn in fcrn])

    log.debug('Testing:')
    log.debug('Original formal CRN:')
    [log.debug('  {}'.format(r)) for r in pretty_crn(fcrn)]
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
    intrp = {deformalize(k, formals): v for k, v in interpretation.items()
                                                } if interpretation else {}
    log.debug('Internal interpretation:')
    [log.debug('  {}: {}'.format(k,v)) for (k, v) in intrp.items()]

    if permissive_depth:
        permissive = [permissive, permissive_depth]
    log.debug(f'Permissive argument: {permissive}')

    out = find_bisimulation(fcrn, icrn, 
                            set(formals), 
                            intrp, 
                            0, 
                            permissive,
                            [{}, 0, [[Counter(),Counter()],Counter()]])

    correct = next(out)
    if correct:
        for intrpout in out:
            yield formalize(intrpout)
    else:
        wintr, max_depth, permissive_failure = next(out)
        wintr = formalize(wintr)
        #max_depth: if > 0, search depth in Qing's algorithm at which intrp was found
        #         if -1, permissive condition was proven false for intrp
        #         if -2, permissive condition could not be proven true for intrp
        #         with specified permissive_depth
        #permissive_failure: if max_depth < 0, permissive_failure[0] is formal
        #    reaction for which permissive condition failed if so and method was
        #    not 'graphsearch', permissive_failure[1] is implementation state
        #    which could not implement the reaction 
        if max_depth >= 0:
            log.info("Delimiting condition cannot be satisfied:")
            if max_depth >= len(fcrn):
                log.info("An implementation reaction cannot be found in the formal CRN.")
            else:
                log.info("A formal reaction cannot be found in the implementation CRN.")
            log.info(f"Maximum search depth reached: {max_depth}")
            log.debug(f"Returning wrong interpretation:\n {wintr}")
        else:
            [[R, P], istate] = permissive_failure
            log.info("Permissive condition cannot be satisfied:")
            log.info("The formal reaction: {}".format(
                '{} -> {}'.format(' + '.join(R), ' + '.join(P))))
            log.info(f"is not possible from implementation state: {istate}")
            if max_depth == -2:
                log.info(f"The maximum trivial reaction chain length {permissive_depth} has been reached.")
        yield None

def modular_crn_bisimulation_test(fcrns, icrns, formals, 
                                  interpretation = None, 
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
        for bisim in crn_bisimulations(fcrn, icrn, 
                                       interpretation = minter, 
                                       formals = mfs, 
                                       permissive = permissive, 
                                       permissive_depth = permissive_depth):
            if bisim is None:
                break
            # Get all formal and implementation species that are in
            # common with at least one other module.
            fsc = {f for f, m in fspc.items() if e in m and len(m) > 1}
            isc = {i for i, m in ispc.items() if e in m and len(m) > 1}
            if is_modular(bisim, icrn, isc, fsc):
                #TODO: re-insert the iterate part here ...
                inter.update(bisim)
                break
            log.debug(f'Skipping non-modular bisimulation: {bisim}')
        else:
            return [False, None]
    return [True, inter]

def crn_bisimulation_test(fcrn, icrn, formals, 
                          interpretation = None,
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
                                 permissive = permissive,
                                 permissive_depth = permissive_depth)

    bisim = next(iterator)
    if bisim is not None:
        return [True, bisim]
    else:
        return [False, None]
