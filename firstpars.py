#!/usr/bin/env python2

import sys
from  igraph import *
import re
g = Graph()
f = open(sys.argv[1], 'r')
w = open(sys.argv[1]+".legend" , 'w')
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
g.add_vertices(nodecount+1)
edgecount = 0
while True:
    raw = f.readline().strip()
    if not line: break;
    if raw == '@attributes':
        break
    el = raw.split()
    g.add_edges([(int(el[0]),int(el[1]))])
    if  el[5]=="4":
        g.es[edgecount]["type"] = 4
        g.es[edgecount]["label"] = el[2]
    elif el[5]=="2":
        g.es[edgecount]["type"] = 2
        g.es[edgecount]["label"] = el[2]
        g.es[edgecount]["Exon"] = el[3][1:]
        # todo does this work???
    else:
        #w.write(str(i)+"\t"+str(el[3][1:])+"\n")
        #print el[0],"->",el[1], '[label="', "S"+str(i), el[6] ,'"];';
        g.es[edgecount]["type"] = 1
        g.es[edgecount]["label"] = el[2]
        g.es[edgecount]["Flow/Capacity"] = el[6]
        g.es[edgecount]["binExon"] = el[4][:-1]
    edgecount = edgecount + 1
print g.get_edgelist()
