__author__ = 'priya'
import numpy as np
import math
import random
# import time
import csv
import sys
import numbers
import collections

## First pass ##
def firstpass_calcdegree(graphpath, thisdelimiter):
    degrees={}

    reader = csv.reader(open(graphpath, 'rb'), delimiter=thisdelimiter)
    for x in reader:
        u = x[0]; v = x[1]

        if degrees.has_key(u):
            degrees[u]=degrees[u]+1
        else:
            degrees[u]=1

        if degrees.has_key(v):
            degrees[v]=degrees[v]+1
        else:
            degrees[v]=1

    spaceused = getsize(degrees)
    return [degrees, spaceused]


## Second and third passes ##
def secondthirdpass_samplewedges_counttraingles_observedtriangles_FAST(graphpath, degrees,thisdelimiter, true_corenumbers):
    #  "secondthirdpass_samplewedges_counttraingles"
    low_degree=10
    num_sampled_ngb={};
    num_sampled_ngb_new={}; lowerbound_corenumbers={}
    for node in degrees.keys():
        num_sampled_ngb[node] = min(int(np.ceil(float(degrees[node])/2)), 50)
        num_sampled_ngb_new[node]=0
    ngbwatchlist={}

    reader = csv.reader(open(graphpath, 'rb'), delimiter=thisdelimiter)
    for x in reader:
        u = x[0]; v = x[1]

        if degrees[u] > low_degree:
            if random.random() <= float(num_sampled_ngb[u])/degrees[u]:
                if v in ngbwatchlist:
                    ngbwatchlist[v].append(u)
                else:
                    ngbwatchlist[v]=[u]
                num_sampled_ngb_new[u]+=1

        if degrees[v] > low_degree:
            if random.random() <= float(num_sampled_ngb[v])/degrees[v]:
                if u in ngbwatchlist:
                    ngbwatchlist[u].append(v)
                else:
                    ngbwatchlist[u]=[v]
                num_sampled_ngb_new[v]+=1

    tricount={}
    for node in degrees.keys():
        num_sampled_ngb[node] = num_sampled_ngb_new[node]
        tricount[node] = 0

    reader = csv.reader(open(graphpath, 'rb'), delimiter=thisdelimiter)
    for x in reader:
        u = x[0]; v = x[1]
        if u not in ngbwatchlist or v not in ngbwatchlist :
            continue

        ngbu = ngbwatchlist[u]
        ngbv = ngbwatchlist[v]
        int_set = list(set(ngbu).intersection(set(ngbv)))
        for node in int_set:
            tricount[node] = tricount[node] + 1

    for node in tricount.keys():
        if degrees[node] <= low_degree:
            lowerbound_corenumbers[node]=1
        else:
            observed_triangles = tricount[node]
            lowerbound_triangles = observed_triangles
            edges_in_egonet = lowerbound_triangles+degrees[node]
            minmaxcliquesize = max_size_of_largest_clique(degrees[node]+1, edges_in_egonet)
            lowerbound_corenumbers[node] = min(max((minmaxcliquesize-1), 1), degrees[node])

    return lowerbound_corenumbers


## A NimbleCore pass ##
def asinglepass(graphpath, degrees2, corenumbers_old, thisdelimiter, lowerbound_estimates_old, binvect_lookup, binvect_lookup_lower):
    bincounts_dict={}; binvals_dict={}
    bincounts_dict_lower={}; binvals_dict_lower={}

    ## Initializing binvals and bincounts ##
    for node in degrees2.keys():
        last_corenumber=corenumbers_old[node]
        binvals_dict[node] = create_binsvals_for_for_a_node(last_corenumber)
        bincounts_dict[node] = np.zeros(len(binvals_dict[node]))

        binvals_dict_lower[node] = create_binsvals_for_for_a_node(last_corenumber)
        bincounts_dict_lower[node] = np.zeros(len(binvals_dict_lower[node]))

    ## A pass ##
    reader = csv.reader(open(graphpath, 'rb'), delimiter=thisdelimiter)
    for x in reader:
        u = x[0]; v = x[1]

        ################################################
        thisnode = u; ngb = v;
        if degrees2[thisnode]>1:
            corenumbers_old_thisnode = corenumbers_old[thisnode]
            corenumbers_old_ngb =corenumbers_old[ngb]
            bincounts_dict, binvals_dict = update_bins_given_node(thisnode, ngb, bincounts_dict, binvals_dict, corenumbers_old_ngb, corenumbers_old_thisnode, binvect_lookup )

        ################################################
        thisnode = v; ngb = u;
        if degrees2[thisnode]>1:
            # print degrees2[ngb]
            corenumbers_old_thisnode = corenumbers_old[thisnode]
            corenumbers_old_ngb =corenumbers_old[ngb]
            bincounts_dict, binvals_dict = update_bins_given_node(thisnode, ngb, bincounts_dict, binvals_dict, corenumbers_old_ngb, corenumbers_old_thisnode, binvect_lookup )

        ################################################
        ################################################

        thisnode = u; ngb = v;
        if degrees2[thisnode]>1:
            lowerbound_estimates_old_thisnode = lowerbound_estimates_old[thisnode]
            lowerbound_estimates_old_ngb =lowerbound_estimates_old[ngb]
            corenumbers_old_thisnode = corenumbers_old[thisnode]
            bincounts_dict_lower, binvals_dict_lower = update_bins_LOWERBOUND_given_nodes(thisnode, ngb, bincounts_dict_lower, binvals_dict_lower, lowerbound_estimates_old_ngb, lowerbound_estimates_old_thisnode, corenumbers_old_thisnode, binvect_lookup_lower)

        ################################################
        thisnode = v; ngb = u;
        if degrees2[thisnode]>1:
            lowerbound_estimates_old_thisnode = lowerbound_estimates_old[thisnode]
            lowerbound_estimates_old_ngb =lowerbound_estimates_old[ngb]
            corenumbers_old_thisnode = corenumbers_old[thisnode]
            bincounts_dict_lower, binvals_dict_lower = update_bins_LOWERBOUND_given_nodes(thisnode, ngb, bincounts_dict_lower, binvals_dict_lower, lowerbound_estimates_old_ngb, lowerbound_estimates_old_thisnode, corenumbers_old_thisnode, binvect_lookup_lower)

    corenumbers_new={}; lowerbound_estimates_new={}

    ## Calculating new estimates from new binvals and bincounts ##
    for node in degrees2.keys():
        if degrees2[node]==1:
            corenumbers_new[node] =1
            lowerbound_estimates_new[node] = 1
        else:
            last_corenumber= corenumbers_old[node]
            corenumbers_new[node] = getHindex_from_bincounts(bincounts_dict[node],corenumbers_old[node], binvals_dict[node])

            last_corenumber_lowerbound= lowerbound_estimates_old[node]
            lowerbound_estimates_new[node] = getHindex_from_bincounts_lowerbound(bincounts_dict_lower[node],lowerbound_estimates_old[node], binvals_dict_lower[node], corenumbers_new[node] )

    spaceused=getsize(bincounts_dict)+getsize(binvals_dict)+\
              getsize(bincounts_dict_lower)+getsize(binvals_dict_lower)+\
              getsize(corenumbers_old)+getsize(lowerbound_estimates_old)

    return corenumbers_new, spaceused, lowerbound_estimates_new
#######################


## lower bound of the max clique in the egonet of a node ##
def max_size_of_largest_clique(num_nodes, num_edges):
    r=1;
    limit_upper = num_nodes-1;
    while num_edges >= limit_upper and r <= num_nodes:
        r=r+1
        limit_upper = np.floor((r-1)*(num_nodes*num_nodes)/float(r*2))
    return r-1

#######################################################
################# helper functions ########################

## Create binvals ##
def create_binsvals_for_for_a_node(last_corenumber):
    numbins = int(math.ceil(np.log2(last_corenumber))+1)
    normal_binvals = np.zeros(numbins)
    for i in range(0,numbins-1):
        normal_binvals[i] = pow(2,i)
    normal_binvals[numbins-1] = last_corenumber

    NC_binvals = np.zeros(numbins)
    NC_binvals[0] = 1
    NC_binvals[-1] = last_corenumber

    for i in range(1,numbins-1):
        binid=numbins-1-i
        # print binid
        diff = normal_binvals[i] - normal_binvals[i-1]
        NC_binvals[binid] = NC_binvals[binid+1] - diff

    return NC_binvals


## Calculate which bin to add a value to ##
def calc_bin_to_add_to(val, ngb_val):
    if ngb_val==1:
        return 0
    if ngb_val>val:
        return int(math.ceil(np.log2(val))+1)-1
    else:
        numbins = int(math.ceil(np.log2(val))+1)
        intans = math.floor(np.log2(val-ngb_val+1))
        ans = numbins-1 - intans
        return ans

## Calculate which bin to add a value to; for lower bound estimate ##
def calc_bin_to_add_to_lower(val, ngb_val):
    if ngb_val==1:
        return 0
    if ngb_val>val:
        return int(math.ceil(np.log2(val))+1)-1
    else:
        numbins = int(math.ceil(np.log2(val))+1)
        intans = math.ceil(np.log2(val-ngb_val+1))
        ans = numbins-1 - intans
        return ans



def update_bins_given_node(thisnode, ngb, bincounts_dict, binvals_dict, corenumbers_old_ngb, corenumbers_old_thisnode, binvect_lookup):
    if corenumbers_old_ngb >= corenumbers_old_thisnode:
        bincounts_dict[thisnode][-1]= bincounts_dict[thisnode][-1]+1
    else:
        bin_to_add_to = binvect_lookup[corenumbers_old_thisnode][corenumbers_old_ngb]
        bincounts_dict[thisnode][bin_to_add_to] +=1
    return bincounts_dict, binvals_dict

def update_bins_LOWERBOUND_given_nodes(thisnode, ngb, bincounts_dict, binvals_dict, lowerbound_estimates_old_ngb, lowerbound_estimates_old_thisnode, corenumbers_old_thisnode, binvect_lookup_lower):
    if lowerbound_estimates_old_ngb >= corenumbers_old_thisnode:
        bincounts_dict[thisnode][-1]= bincounts_dict[thisnode][-1]+1
    else:
        bin_to_add_to = binvect_lookup_lower[corenumbers_old_thisnode][lowerbound_estimates_old_ngb]
        bincounts_dict[thisnode][bin_to_add_to] +=1
    return bincounts_dict, binvals_dict


## Calculate estimated error from the upper and lower bounds of the core number estimates ##
def estimateerror(corenumbers_estimates,lowerbound):
    esterrlist=[]
    for node in corenumbers_estimates.keys():
        estimated_error = abs(float(corenumbers_estimates[node])-float(lowerbound[node]))/float(lowerbound[node])
        esterrlist.append(estimated_error)
    esterr = np.mean(esterrlist)
    return esterr


def getHindex_from_bincounts(bincounts, last_corenumber, binvals):
    numbins=len(bincounts)
    idx=0; maxofmins=0
    node_degree= sum(bincounts)
    for c in range(numbins):
        if bincounts[c]>0:

            d_minus_i_plus_1_end= node_degree - (idx+1) + 1 # 2 - 2 +2 = 2
            idx+= bincounts[c]
            d_minus_i_plus_1_start = node_degree - idx + 1 # 2 - 2 +1 = 1

            for r in range(int(d_minus_i_plus_1_start), int(d_minus_i_plus_1_end+1)):
                d_minus_i_plus_1 = r
                if maxofmins < min(d_minus_i_plus_1, binvals[c]):
                    maxofmins = min(d_minus_i_plus_1, binvals[c])
    maxofmins = min(maxofmins, last_corenumber)
    return maxofmins

def getHindex_from_bincounts_lowerbound(bincounts, last_lowebound, binvals, last_corenumber):

    if last_lowebound > last_corenumber:
        last_lowebound = last_corenumber

    numbins=len(bincounts)
    idx=0; maxofmins=0
    node_degree= sum(bincounts)

    for c in range(numbins):
        if bincounts[c]>0:
            d_minus_i_plus_1_end= node_degree - (idx+1) + 1 # 2 - 2 +2 = 2
            idx+= bincounts[c]
            d_minus_i_plus_1_start = node_degree - idx + 1 # 2 - 2 +1 = 1
            for r in range(int(d_minus_i_plus_1_start), int(d_minus_i_plus_1_end+1)):
                d_minus_i_plus_1 = r
                if maxofmins < min(d_minus_i_plus_1, binvals[c]):
                    maxofmins = min(d_minus_i_plus_1, binvals[c])

    if maxofmins > last_corenumber:
        maxofmins = last_corenumber
    else:
        maxofmins = max(maxofmins, last_lowebound)
    return maxofmins


## For efficient lookup ##
def build_binvalue_hashtable_fromcalc( degrees):
    binvect_lookup = {}
    binvect_lookup_lower = {}

    for val in list(set(degrees.values())):
        binvect_lookup[val]={}
        binvect_lookup_lower[val]={}
        for ngb_val in range(1, val+1):
            binvect_lookup[val][ngb_val] = calc_bin_to_add_to(val, ngb_val)
            binvect_lookup_lower[val][ngb_val] = calc_bin_to_add_to_lower(val, ngb_val)

    return binvect_lookup, binvect_lookup_lower



#######################################################
#######################################################


## For efficient lookup ##
def update_binvalue_hashtable_fromcalc(corenumbers_estimates, lowerbound_estimates, binvect_lookup, binvect_lookup_lower):
    max_cnest = max(corenumbers_estimates.values())
    max_cnest_low = max(lowerbound_estimates.values())
    for key in binvect_lookup_lower.keys():
        if key> max_cnest_low:
            binvect_lookup_lower.pop(key)
    for key in binvect_lookup.keys():
        if key> max_cnest:
            binvect_lookup.pop(key)

    for val in list(set(corenumbers_estimates.values())):
        val = int(val)
        if val not in corenumbers_estimates:
            binvect_lookup[val]={}
            for ngb_val in range(1, val+1):
                binvect_lookup[val][ngb_val] = calc_bin_to_add_to(val, ngb_val)

    for val in list(set(corenumbers_estimates.values())):
        val = int(val)
        if val not in lowerbound_estimates:
            binvect_lookup_lower[val]={}
            for ngb_val in range(1, val+1):
                binvect_lookup_lower[val][ngb_val] = calc_bin_to_add_to_lower(val, ngb_val)

    return binvect_lookup, binvect_lookup_lower


## Calculate the error of the estimates ##
def evaluate_cnestimates(corenumbers_true, corenumbers_estimates, degrees):
    high_degree=math.ceil(np.percentile(degrees.values(),90))
    error_high=[]; errors=[]

    for node in corenumbers_estimates.keys():
        thiserr=abs(float(corenumbers_estimates[node])-float(corenumbers_true[node]))/float(corenumbers_true[node])
        errors.append(thiserr)
        if degrees[node]>high_degree:
            error_high.append(thiserr)

    error = np.mean(errors)
    if not error_high :
        err_high=0
    else:
        err_high= np.mean(error_high)

    return [error, err_high]


## Get size of objects ##
# From http://stackoverflow.com/questions/449560/how-do-i-determine-the-size-of-an-object-in-python
def getsize(obj):
    # recursive function to dig out sizes of member objects:
    def inner(obj, _seen_ids = set()):
        obj_id = id(obj)
        if obj_id in _seen_ids:
            return 0
        _seen_ids.add(obj_id)
        size = sys.getsizeof(obj)
        if isinstance(obj, (basestring, numbers.Number, xrange)):
            pass # bypass remaining control flow and return
        elif isinstance(obj, (tuple, list, set, frozenset)):
            size += sum(inner(i) for i in obj)
        elif isinstance(obj, collections.Mapping) or hasattr(obj, 'iteritems'):
            size += sum(inner(k) + inner(v) for k, v in obj.iteritems())
        else:
            attr = getattr(obj, '__dict__', None)
            if attr is not None:
                size += inner(attr)
        return size
    return inner(obj)