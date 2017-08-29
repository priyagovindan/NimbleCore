'''
Created on Dec 17, 2014

@author: priya
'''

import numpy as np
import kcore_decomposition as kdecomp
import csv
import gc
import NimbleCore_functions as nc
import timeit
import networkx as nx
import sys

def _template_func(setup, func):
    """Create a timer function. Used if the "statement" is a callable."""
    def inner(_it, _timer, _func=func):
        setup()
        _t0 = _timer()
        for _i in _it:
            retval = _func()
        _t1 = _timer()
        return _t1 - _t0, retval
    return inner
timeit._template_func = _template_func


def main(graphpath, thisdelimiter):

    print
    print 'Input graph is ',graphpath

    G=nx.read_edgelist(graphpath,delimiter=thisdelimiter, nodetype=str)

    ## In-memory k-core decomposition ##
    true_corenumbers = kdecomp.getCoreNumbers(G)

    print 'Number of nodes = ', G.number_of_nodes()
    print 'Number of edges = ', G.number_of_edges()

    kcoredecomp_space = float(nc.getsize(G))/(1024*1024)
    del G

    ## NimBleCore ##
    totaltime = 0; allspaceused = []

    ## First pass: Calculate degree, the first upper bound estimates ##
    t1,[degrees, spaceused ]= timeit.Timer(lambda: nc.firstpass_calcdegree(graphpath, thisdelimiter)).timeit(number=1)

    ## Computing bin lookups for speed efficiency. This can be pre-computed and stored in a file, for further speedup ##
    binvect_lookup, binvect_lookup_lower = nc.build_binvalue_hashtable_fromcalc(degrees)

    ## Second and third passes: Compute lower bound estimates for the nodes ##
    t2,lowerbound_estimates= timeit.Timer(lambda: nc.secondthirdpass_samplewedges_counttraingles_observedtriangles_FAST(graphpath, degrees,thisdelimiter, true_corenumbers)).timeit(number=1)

    totaltime += t1+t2

    ## Estimate error based on estimates from first 3 passes ##
    esterr= nc.estimateerror(degrees, lowerbound_estimates)

    error_estq_OLD = esterr

    corenumbers_old = degrees.copy()
    lowerbound_estimates_old = lowerbound_estimates.copy()

    totalpasses=3; flag=1

    while flag:

        totalpasses= totalpasses+ 1

        ## A NimbleCore pass ##
        timetaken_1pass,[corenumbers_estimates, spaceused, lowerbound_estimates]=\
            timeit.Timer(lambda: nc.asinglepass(graphpath, degrees, corenumbers_old.copy(), thisdelimiter, lowerbound_estimates_old.copy(), binvect_lookup, binvect_lookup_lower)).timeit(number=1)

        allspaceused.append(spaceused)
        totaltime += timetaken_1pass

        ## Estimate error ##
        esterr= nc.estimateerror(corenumbers_estimates, lowerbound_estimates)

        corenumbers_old = corenumbers_estimates.copy()
        lowerbound_estimates_old = lowerbound_estimates.copy()

        ## Checking the stopping condition ##
        if error_estq_OLD - esterr < 0.01:
            flag = 0
        error_estq_OLD = esterr

        ## Computing bin lookups for speed efficiency. This can be pre-computed and stored in a file, for further speedup ##
        binvect_lookup, binvect_lookup_lower = nc.update_binvalue_hashtable_fromcalc(corenumbers_estimates, lowerbound_estimates, binvect_lookup, binvect_lookup_lower)

        gc.collect()

    ## End of all passes ##

    [error, err_high] = nc.evaluate_cnestimates(true_corenumbers, corenumbers_estimates, degrees)

    spaceused = float(max(allspaceused))/(1024*1024)

    ## Writing core number estimates to file ##
    writer = csv.writer(open('NimbleCore_corenumbers_'+graphpath,'wb'))
    writer.writerow(['node', 'corenumber'])
    for key, value in corenumbers_estimates.items():
        writer.writerow([key, value])

    ## Writing true core numbers to file #
    writer = csv.writer(open('true_corenumbers_'+graphpath,'wb'))
    writer.writerow(['node', 'corenumber'])
    for key, value in true_corenumbers.items():
        writer.writerow([key, value])

    print
    print 'Relative Error = ', error
    print 'Relative Error of top 10% high degree nodes = ', err_high

    print
    print 'Space used (in megabytes) = ',spaceused
    print 'Space used by in-memory k-core decomposition (in megabytes) = ',kcoredecomp_space

    print
    print 'Number of passes = ', totalpasses
    print 'Time taken (in seconds) = ', totaltime

    print
    print 'The true and estimated core numbers have been saved in the current directory.'

    ##


    gc.collect()


##### END #####



graphpath = sys.argv[1]

thisdelimiter= ' '

main(graphpath, thisdelimiter)

