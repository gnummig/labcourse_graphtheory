#!/usr/bin/env python2

import sys
from  igraph import *
import re

def checkDegree(x):
	if x<2:
		return 0
	else:
		return 1

#parse original graph
def getGraph():
    graph = Graph(directed=True)
    while True:
        line = f.readline().strip()
        if not line:
            break;
        if line == "@nodes":
            f.readline()
            break;
    while True:
        line = f.readline().strip()
        if not line: break;
        if line == "@arcs":
            f.readline()
            break;
        nodecount=int(line)
    graph.add_vertices(nodecount+1)
    edgecount = 0
    while True:
        raw = f.readline().strip()
        if not line: break;
        if raw == '@attributes':
            break;
        el = raw.split()
        graph.add_edges([(int(el[0]),int(el[1]))])
        if  el[5]=="4":
            graph.es[edgecount]["type"] = 4
            graph.es[edgecount]["label"] = el[2]
        if el[5]=="2":
            graph.es[edgecount]["type"] = 2
            graph.es[edgecount]["label"] = el[2]
            graph.es[edgecount]["Exon"] = el[3][1:]
            # todo does this work???
        else:
            #w.write(str(i)+"\t"+str(el[3][1:])+"\n")
            #print el[0],"->",el[1], '[label="', "S"+str(i), el[6] ,'"];';
            graph.es[edgecount]["type"] = 1
            graph.es[edgecount]["label"] = el[2]
            graph.es[edgecount]["Flow/Capacity"] = el[6]
            graph.es[edgecount]["binExon"] = el[4][:-1]
        edgecount = edgecount + 1
    return graph

f = open(sys.argv[1], 'r')
origGraph = getGraph()
compGraph = getGraph()
resGraph = getGraph()
#w = open(sys.argv[1]+".legend" , 'w')

print "Edgelist of original graph"
print origGraph.get_edgelist()
print "Edgelist of compacted graph"
print compGraph.get_edgelist()
print "Edgelist of resolved graph"
print resGraph.get_edgelist()

problemIn = map(checkDegree, compGraph.degree(type="in"))
problemOut = map(checkDegree, compGraph.degree(type="out"))
problemNodes = [a*b for a,b in zip(problemIn, problemOut)]
print "Problem Nodes in compacted graph:"
print problemNodes
