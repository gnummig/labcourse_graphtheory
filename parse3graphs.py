#!/usr/bin/env python2

import sys
from  igraph import *
import re

def checkDegree(x):
	if x<2:
		return 0
	else:
		return 1

origGraph = Graph(directed=True)
compGraph = Graph(directed=True)
resGraph = Graph(directed=True)
f = open(sys.argv[1], 'r')
#w = open(sys.argv[1]+".legend" , 'w')

#parse original graph
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
origGraph.add_vertices(nodecount+1)
edgecount = 0
while True:
    raw = f.readline().strip()
    if not line: break;
    if raw == '@attributes':
        break
    el = raw.split()
    origGraph.add_edges([(int(el[0]),int(el[1]))])
    if  el[5]=="4":
        origGraph.es[edgecount]["type"] = 4
        origGraph.es[edgecount]["label"] = el[2]
    elif el[5]=="2":
        origGraph.es[edgecount]["type"] = 2
        origGraph.es[edgecount]["label"] = el[2]
        origGraph.es[edgecount]["Exon"] = el[3][1:]
        # todo does this work???
    else:
        #w.write(str(i)+"\t"+str(el[3][1:])+"\n")
        #print el[0],"->",el[1], '[label="', "S"+str(i), el[6] ,'"];';
        origGraph.es[edgecount]["type"] = 1
        origGraph.es[edgecount]["label"] = el[2]
        origGraph.es[edgecount]["Flow/Capacity"] = el[6]
        origGraph.es[edgecount]["binExon"] = el[4][:-1]
    edgecount = edgecount + 1
print "Edgelist of original graph"
print origGraph.get_edgelist()

#parse compacted graph
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
compGraph.add_vertices(nodecount+1)
edgecount = 0
while True:
    raw = f.readline().strip()
    if not line: break;
    if raw == '@attributes':
        break
    el = raw.split()
    compGraph.add_edges([(int(el[0]),int(el[1]))])
    if  el[5]=="4":
        compGraph.es[edgecount]["type"] = 4
        compGraph.es[edgecount]["label"] = el[2]
    elif el[5]=="2":
        compGraph.es[edgecount]["type"] = 2
        compGraph.es[edgecount]["label"] = el[2]
        compGraph.es[edgecount]["Exon"] = el[3][1:]
        # todo does this work???
    else:
        #w.write(str(i)+"\t"+str(el[3][1:])+"\n")
        #print el[0],"->",el[1], '[label="', "S"+str(i), el[6] ,'"];';
        compGraph.es[edgecount]["type"] = 1
        compGraph.es[edgecount]["label"] = el[2]
        compGraph.es[edgecount]["Flow/Capacity"] = el[6]
        compGraph.es[edgecount]["binExon"] = el[4][:-1]
    edgecount = edgecount + 1
print "Edgelist of compacted graph"
print compGraph.get_edgelist()

#parse resolved graph
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
resGraph.add_vertices(nodecount+1)
edgecount = 0
while True:
    raw = f.readline().strip()
    if not line: break;
    if raw == '@attributes':
        break
    el = raw.split()
    resGraph.add_edges([(int(el[0]),int(el[1]))])
    if  el[5]=="4":
        resGraph.es[edgecount]["type"] = 4
        resGraphes[edgecount]["label"] = el[2]
    elif el[5]=="2":
        resGraph.es[edgecount]["type"] = 2
        resGraph.es[edgecount]["label"] = el[2]
        resGraph.es[edgecount]["Exon"] = el[3][1:]
        # todo does this work???
    else:
        #w.write(str(i)+"\t"+str(el[3][1:])+"\n")
        #print el[0],"->",el[1], '[label="', "S"+str(i), el[6] ,'"];';
        resGraph.es[edgecount]["type"] = 1
        resGraph.es[edgecount]["label"] = el[2]
        resGraph.es[edgecount]["Flow/Capacity"] = el[6]
        resGraph.es[edgecount]["binExon"] = el[4][:-1]
    edgecount = edgecount + 1
print "Edgelist of resolved graph"
print resGraph.get_edgelist()

problemIn = map(checkDegree, compGraph.degree(type="in"))
problemOut = map(checkDegree, compGraph.degree(type="out"))
problemNodes = [a*b for a,b in zip(problemIn, problemOut)]
print "Problem Nodes in compacted graph:"
print problemNodes



