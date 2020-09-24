from collections import Counter

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
    from .crn_bisimulation import msleq, interpret

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


