from collections import Counter
from itertools import product, permutations

# helper functions copied from crn_bisimulation.py
def msleq(x,y):
    # True if (multisets) x <= y (vector comparison)
    for k in x:
        if x[k] > y[k]:
            return False
    return True

def mstimes(s, l):
    # return l*s for integer l, multiset s
    c = Counter()
    for k in s:
        c[k] = l * s[k]
    return c

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



# DEPRECATED STATS HERE
def find_one_trivial(fcrn, icrn, fs, intrp, permcheck = 'graphsearch'):

    # All unknown implementation reactions (i.e. those with some unknown
    # species) must be trivial.  Build the matrix for the the "solve" function,
    # to see whether the interpretation can be completed as required.  Also
    # note that "unknow" is not used; the unknown species are recalculated here
    # because "unknow" contains only those implementation reactions that have
    # not been solved by the row search, but other rows (i.e. implementation
    # reactions) may be implicitly solved by the partial interpretation.
    sicrn = subst(icrn, intrp)

    # List of implementation reactions with an unknown species.
    unknown = [i for i, irxn in enumerate(sicrn) if len(set(irxn[0]) - fs) or \
                                                   len(set(irxn[1]) - fs)]
    # All species that remain unknown in current (partial) interpretation.
    unassigned = set(s for irxn in sicrn for s in set(irxn[0]) - fs | set(irxn[1]) - fs)
    
    # Find species that do not satisfy atomic condition
    atoms_needed = list(fs - set().union(*[set(a) for a in intrp.values() if sum(a.values()) == 1]))

    log.info(f'Trivial reaction solver: {unknown=}')
    log.info(f'Missing atoms: {atoms_needed}')
    log.info(f'Unknown species: {unassigned}')

    for assign in permutations(unassigned, len(atoms_needed)): # works even if l == 0
        # each assign is a tuple of implementation species to be interpreted as
        # exactly one formal species, matching the order of atoms_needed.
        log.debug(f'{assign=}')

        parti = intrp.copy()
        for i, a in enumerate(atoms_needed):
            assert assign[i] not in parti
            parti[assign[i]] = Counter({a: 1})
        log.debug(f'{parti=}')

        sicrn = subst(icrn, parti)
        T = updateT(fcrn, sicrn, fs)
        if (not checkT(T)) or any([T[i][-1] is False for i in unknown]):
            continue

        ulist = [sp for sp in unassigned if sp not in assign]
        log.debug(f'{ulist=}')
        if not ulist:
            out, info = check_permissive(fcrn, icrn, fs, parti, permcheck)
            if out:
                yield parti
            continue

        # A mini table. 
        # For every unknown implementation reaction make a list:
        #   - for each unassigend species store #ur - #up
        a = []
        for i in unknown:
            irxn = sicrn[i]
            log.debug(irxn)
            # A list of #ur - #up
            a.append([irxn[0][u]-irxn[1][u] for u in ulist] + [0])
        log.debug(f'{a=}')

        # prepare for adding species to the counter ...
        for u in ulist:
            parti[u] = Counter()
        log.debug(f'{parti=}')

        check = True
        for fsp in fs:
            log.warning(f'{fsp=}')
            for e, i in enumerate(unknown):
                irxn = sicrn[i]
                # A list of #fsr - #fsp
                a[e][-1] = irxn[0][fsp]-irxn[1][fsp]
            log.debug(f'{a=}')

            # Alright, this can find one, but only one solution... 
            s = solve(a)
            log.debug(f'{s=}')
            if s == []:
                check = False
                break
            else:
                for j, u in enumerate(ulist):
                    parti[u][fsp] = s[j]
        if check:
            for isp in ulist:
                parti[isp] = parti[isp] + Counter()
            out, info = check_permissive(fcrn, icrn, fs, parti, permcheck)
            if out:
                yield parti
    return

def moduleCond(module, formCommon, implCommon, intrp):
    """
    # check whether the modularity condition (every implementation species can
    # turn into common species with the same interpretation) is satisfied
    # assumes intrp is complete and filtered, so intrp.keys() is a list of all
    # implementation species in this module (and no others) algorithm is
    # basically the graphsearch algorithm from perm, where all "minimal states"
    # are each exactly one implementation species

    # canBreak[k] is:
    #   True if species k is known to decompose (via trivial reactions) into
    #     species which are each either in implCommon or have an interpretation
    #     containing nothing in formCommon, or known to decompose as such if
    #     some set of other species each with interpretation strictly < do
    #   [reach,produce] where reach is set of implementation species with
    #     interpretation equal to that of k (known) reachable from k, and
    #     produce is set of null species producible in a loop from k to k
    """
    module = [[Counter(part) for part in rxn] for rxn in module]
    intrp = {k : Counter(v) for k, v in intrp.items()}

    canBreak = {k: ((k in implCommon) or set(intrp[k]).isdisjoint(formCommon)
                    or [set(),set()])
                for k in intrp}
    changed = True

    tr = [rxn for rxn in module if (lambda x,y: msleq(x,y) and msleq(y,x))
          (interpret(rxn[0],intrp),interpret(rxn[1],intrp))]

    while changed:
        changed = False
        for k in canBreak:
            if canBreak[k] is True: continue

            for rxn in tr:
                if k in rxn[0] and \
                   set(rxn[0] - Counter([k])).issubset(canBreak[k][1]):
                    # reactants of rxn are one copy of k and some null species
                    #  producible in a loop from k to k
                    nulls = set()
                    theOne = None
                    for sp in rxn[1]:
                        if intrp[sp] == Counter():
                            nulls.add(sp)
                        elif not msleq(intrp[k],intrp[sp]):
                            canBreak[k] = True
                            changed = True
                            break
                        else:
                            if canBreak[sp] is True:
                                canBreak[k] = True
                                changed = True
                                break
                            theOne = sp

                    if canBreak[k] is True:
                        break

                    if theOne not in canBreak[k][0]:
                        canBreak[k][0].add(theOne)
                        changed = True

                    if not (canBreak[theOne][0] <= canBreak[k][0]):
                        canBreak[k][0] |= canBreak[theOne][0]
                        changed = True

                    if k in canBreak[theOne][0]:
                        loopable = nulls | canBreak[theOne][1]
                        if not loopable <= canBreak[k][1]:
                            canBreak[k][1] |= loopable
                            changed = True

    return all([canBreak[k] is True for k in canBreak])

def update(fcrn, icrn, fs):
    # the um, er, um, completely recalculates the table from scratch.
    # assumes subst has already been applied to implementation icrn.
    # This should be logically equivalent to the UpdateTable in the MS thesis
    # for compiled DNA reactions (we think).

    # WARNING:
    # If an implementation CRN has directly catalytic species, the code below
    # may fail (though thesis psuedocode is correct).  
    # E.g.  i3 + i7 --> i12 + i7 + i8
    m = len(icrn)
    r = []
    for i in range(len(icrn)):
        rr = []
        for j in range(len(fcrn)):
            t1 = fcrn[j][0] - icrn[i][0]
            t2 = icrn[i][0] - fcrn[j][0]
            t3 = fcrn[j][1] - icrn[i][1]
            t4 = icrn[i][1] - fcrn[j][1]
            if set(t2).isdisjoint(set(fs)) and set(t4).isdisjoint(set(fs)):
                if list(t1.keys()) == []:
                    if list(t3.keys()) == []:
                        rr.append(True)
                    else:
                        if list(t4.keys()) == []:
                            rr.append(False)
                        else:
                            rr.append(True)
                else:
                    if list(t2.keys()) == []:
                        rr.append(False)
                    else:
                        if list(t3.keys()) == []:
                            rr.append(True)
                        else:
                            if list(t4.keys()) == []:
                                rr.append(False)
                            else:
                                rr.append(True)
            else:
                rr.append(False)
        t1 = icrn[i][0] - icrn[i][1]
        t2 = icrn[i][1] - icrn[i][0]
        if (set(t1).isdisjoint(set(fs)) or not set(t2).issubset(set(fs))) and \
                (set(t2).isdisjoint(set(fs)) or not set(t1).issubset(set(fs))):
            rr.append(True)
        else:
            rr.append(False)
        r.append(rr)
    return r

def subsets(x):
    # generates all subsets of multiset x
    [ks, vs] = zip(*x.items()) if len(x) else [[],[]]
    vs = [list(range(v + 1)) for v in vs]
    for prod in product(*vs):
        # calls to keys and values with no dictionary modifications in between
        #  should produce the keys and values in the same order,
        #  and product should respect that order (I think)
        yield Counter(dict(list(zip(ks, prod)))) + Counter()

def enum(n, s, weights = None):
    """Partition multiset s into n ordered (possibly empty) parts.  

    For example:
        enum(2, [a, b]) = [ [[], [a, b]], 
                            [[a], [b]], 
                            [[b], [a]], 
                            [[a, b],[]] ]

    If weights are given, enumerates all lists l of n multisets such that 
        s = sum(weights[i] * l[i]) 
    (if weights are not given, equivalent to weights[i] = 1 for all i)
    """
    if weights is None:
        weights = [1] * n

    if n == 0:
        yield []
        return

    def msdiv(s, l):
        # return s/l for integer l, multiset s if l divides s, otherwise False
        # l divides s if l divides s[k] for every key k
        c = Counter()
        for k in s:
            c[k] = int(s[k]/l)
            if c[k] * l != s[k]:
                return False
        return c

    if len(weights) < n:
        raise IndexError(f'{len(weights)} weights given for {n} parts.')
    elif weights[0] < 0:
        raise ValueError('Negative weight given.')
    elif weights[0] == 0:
        for j in enum(n-1, s, weights[1:]):
            yield [Counter()] + j
        return

    if n == 1:
        sdivw = msdiv(s, weights[0])
        if sdivw is not False: 
            yield [sdivw]
        return

    for i in subsets(s):
        ss = mstimes(i, weights[0])
        if not msleq(ss, s):
            continue
        for j in enum(n-1, s-ss):
            yield [i] + j

def searchc(fcrn, icrn, fs, intrp, unknown, depth, permcheck, state):
    # Search column.  I.e. make sure every formal reaction can be implemented.
    """
    fcrn, icrn, fs, unknown and intrp, permcheck remain constant, i think.

    Depth will update, state will update.

    """
    from .crn_bisimulation import subst, updateT, checkT, enumL

    # A quick check if the delimiting condition is still satisfied.
    sicrn = subst(icrn, intrp)
    T = updateT(fcrn, sicrn, fs)
    if not checkT(T):
        yield False
        yield state
        return

    # Rollback to original interpretation if you reach the max_depth,
    # I assume this is a hack to do some sort of depth/breadth first search.
    intr, max_depth = state[0:2]
    if max_depth >= 0 and depth > max_depth:
        intr = intrp.copy()
        max_depth = depth

    found = False # Changes to True if you found something ... but what?

    # Find the next column to solve. 
    c, worst = -1, len(icrn) + 1
    for i in unknown: # The column indices of formal reactions 
        # Count the number of True's in the column.
        tmp = sum(T[j][i] for j in range(len(icrn)))
        if tmp < worst:
            c, worst = i, tmp

    if c < 0: # Done with column search, transition to row search.
        untmp = []
        for i in range(len(icrn)):
            if not (set(sicrn[i][0])-set(fs) == set([]) and \
                    set(sicrn[i][1])-set(fs) == set([])):
                untmp.append(i)
        out = searchr(fcrn, icrn, fs, untmp, intrp, depth, permcheck,
                      [intr, max_depth] + state[2:])
        if next(out):
            yield True
            found = True
            for outintrp in out:
                yield outintrp
        else:
            yield False
            yield next(out)
            return
    else:
        # search_column(T, c, unknown, icrn, fcrn, intrp, sicrn)
        # This part will recursively call searchc 

        untmp = list(unknown)
        untmp.remove(c)
        for k in range(len(icrn)):
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
                    itmp = intrp.copy()
                    itmp.update({k: Counter(v) for k, v in intrpleft.items()})
                    itmp.update({k: Counter(v) for k, v in intrpright.items()})
                    out = searchc(fcrn, icrn, fs, 
                                  itmp, 
                                  untmp, 
                                  depth + 1,
                                  permcheck, 
                                  [intr, max_depth] + state[2:])
                    if next(out):
                        if not found:
                            found = True
                            yield True
                        for outintrp in out:
                            yield outintrp
                    else:
                        state = next(out)
                        intr, max_depth = state[0:2]

    if not found:
        yield False
        yield [intr, max_depth] + state[2:]

def equations(fcrn, icrn, fs, intrp, permcheck, state):
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

def perm(fcrn, icrn, fs, intrp, permcheck, state):
    """ Check the permissive condition.
    
    Uses the original formal CRN  and original implementation CRN
    (without substitutions).

    Args:
        fcrn: The original formal CRN
        icrn: The original implementation CRN
        fs: The formal species.

    """
    from .crn_bisimulation import subst, updateT, checkT, enumL
    if len(permcheck) == 2:
        permissive_depth = permcheck[1]
        permcheck = permcheck[0]
    else:
        permissive_depth = None

    intr, max_depth, permissive_failure = state

    tr = []
    fr = []
    hasht = set([])

    sicrn = subst(icrn, intrp)
    T = updateT(fcrn, sicrn, fs)
    assert checkT(T) # Assuming this has been checked before calling permissive.
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

def permissive(fcrn, icrn, fs, intrp, permcheck):
    """ Check the permissive condition.
    
    Uses the original formal CRN  and original implementation CRN
    (without substitutions).

    Args:
        fcrn: The original formal CRN
        icrn: The original implementation CRN
        fs: The formal species.

    """
    from .crn_bisimulation import subst, updateT, checkT, enumL
    permissive_depth = None

    tr = []
    fr = []
    hasht = set([])

    sicrn = subst(icrn, intrp)
    T = updateT(fcrn, sicrn, fs)
    assert checkT(T) # Assuming this has been checked before calling permissive.
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

def searchr(fcrn, icrn, fs, unknown, intrp, d, permcheck, state, nontriv=False):
    """ Row search: every implementation reaction must interpret to a formal reaction 
        (or be trivial).
    """
    from .crn_bisimulation import subst, updateT, checkT, enumL
    intr, max_depth = state[0:2]
    sicrn = subst(icrn, intrp)
    T = updateT(fcrn, sicrn, fs)
    if not checkT(T):
        # delimiting condition is not satisfied
        yield False
        yield state
        return
    if unknown == []:
        for fsp in fs:
            checkFs = False
            for isp in intrp:
                if intrp[isp] and not intrp[isp] - Counter([fsp]):
                    checkFs = True
                    break
            if not checkFs:
                # atomic condition is not satisfied
                yield False
                yield state
                return
        out = perm(fcrn, icrn, fs, intrp, permcheck, state)
        yield out[0]
        yield out[1]
        return
    if max_depth >= 0 and d > max_depth:
        intr = intrp.copy()
        max_depth = d
    found = False
    min = len(fcrn)+1
    k = -1  # next row reaction we will solve
    nt = 0  # number of possibly trivial reactions according to table
    for i in unknown:
        tmp = T[i].count(True)
        if T[i][-1] == True:
            nt += 1
            tmp -= 1
        if tmp < min and tmp > 0:
            min = tmp
            k = i
    if nt == len(unknown) and not nontriv:  # all unsearched rows could be trivial -- try it!
        out = equations(fcrn, icrn, fs, intrp, permcheck,
                        [intr, max_depth] + state[2:])
        if out[0]:
            yield True
            yield out[1]
            found = True
        else:
            # if we just tried equations and it didn't work, then we
            #  shouldn't try again unless something changes
            nontriv = True
            state = out[1]
            intr, max_depth = state[0:2]
    if k < 0:
        if not found:
            yield False
            yield [intr, max_depth] + state[2:]
        return
    untmp = list(unknown)
    untmp.remove(k)
    if T[k][-1] == True:  # if implementation reaction #k can be trivial, leave it that way and try more rows
        out = searchr(fcrn, icrn, fs, untmp, intrp, d, permcheck,
                      [intr, max_depth] + state[2:], nontriv)
        if next(out):
            if not found:
                found = True
                yield True
            for outintr in out:
                yield outintr
        else:
            state = next(out)
            intr, max_depth = state[0:2]
    n = 0
    for c in range(len(fcrn)): # try to match implementation reaction #k with some formal reaction
        if T[k][c]:        
            ul = sicrn[k][0] - fcrn[c][0]
            kl = list(ul.keys())
            vl = list(ul.values())
            nl = len(kl)
            sl = fcrn[c][0] - sicrn[k][0]
            tmpl = enum(nl, sl, vl)
            ur = sicrn[k][1] - fcrn[c][1]
            kr = list(ur.keys())
            vr = list(ur.values())
            nr = len(kr)
            sr = fcrn[c][1] - sicrn[k][1]
            tmpr = enum(nr, sr, vr)
            for (i,j) in product(tmpl, tmpr):
            # for i in tmpl:
                intrpleft = dict(list(zip(kl, i)))
                # for j in tmpr:
                intrpright = dict(list(zip(kr, j)))
                checkCompatible = True
                for key in intrpleft:
                    if key in intrpright and \
                       any([intrpright[key][fsp] != intrpleft[key][fsp]
                            for fsp in fs]):
                        checkCompatible = False
                        break

                if not checkCompatible:
                    continue

                itmp = intrp.copy()
                itmp.update(intrpleft)
                itmp.update(intrpright)
                out = searchr(fcrn, icrn, fs, untmp, itmp, d+1, permcheck,
                              [intr, max_depth] + state[2:])
                if next(out):
                    if not found:
                        found = True
                        yield True
                    for outintrp in out:
                        yield outintrp
                else:
                    state = next(out)
                    intr, max_depth = state[0:2]
    if not found:
        yield False
        yield [intr, max_depth] + state[2:]
    return

