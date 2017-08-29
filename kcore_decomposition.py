'''
Created on Jan 21, 2015

@author: priya
'''
import networkx as nx

import time



def getCoreNumbers(Gpassed):
    G=Gpassed.copy()
    G=removeSelfLoops(G)
    G=removeSingletons(G)

    degrees=G.degree()
    dmax=max(degrees.values())+2
    corenumbers={}
    
    N=len(degrees.values())
    prevN = N
    for corenumber_k in range(2,dmax):
        t1=time.time()
        whilec=0
        while min(degrees.values())<corenumber_k:
            whilec+=1
            for node in degrees.keys():
                if degrees[node]<corenumber_k:
                    corenumbers[node]=corenumber_k-1
                    G.remove_node(node)

            degrees=G.degree()
            if len(degrees.values())==0 :
                break

        if len(degrees.values())==0 :
            break

    return corenumbers



def removeSelfLoops(G):
    nodes_with_selfloops = G.nodes_with_selfloops()

    for node in nodes_with_selfloops:
        G.remove_edge(node, node)

    return G

def removeSingletons(G):
    degrees=G.degree()

    for node in degrees.keys():
        if degrees[node]==0:
            G.remove_node(node)

    return G

