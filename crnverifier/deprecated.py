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


# DEPRECATED STATS HERE

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
    from .crn_bisimulation import subst, updateT, checkT, enumL, searchr

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

