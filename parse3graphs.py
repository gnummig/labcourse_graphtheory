#!/usr/bin/env python2

import sys
from  igraph import *
import re

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

def printGraphDatatable(graph, name):
    inDegree =  graph.degree(type="in")
    outDegree = graph.degree(type="out")
    problemNodes = [1*(a>1)*(b>1) for a,b in zip(inDegree,outDegree)]
    centrality = graph.as_undirected().evcent()

    #GraphId, VertexCount, Vertex_ID, In_Degree, Out_Degree, ProblemNode?, Centrality
    for i in range(0, graph.vcount()):
        print name +"\t"+ str(graph.vcount()) +"\t" + str(i) + "\t" + str(inDegree[i]) + "\t" + str(outDegree[i]) + "\t" + str(problemNodes[i]) + "\t" + str(centrality[i])
        #print("")
    return

##########
#  main  #
##########


f = open(sys.argv[1], 'r')

# remove first line, after that the exon list followx
f.readline().strip()
while True:
    line = f.readline().strip()
    if not line:
        break;
    if "Multi" in line :
        break;
    exonPos=[[],[]]
    raw = f.readline()
    exonPos[0].append(raw.split( " " )[1].split( "-" )[0])
    exonPos[1].append(raw.split( " " )[1].split( "-" )[1])

origGraph = getGraph()
compGraph = getGraph()
resGraph = getGraph()

if  len(sys.argv)>2:
    printGraphDatatable(origGraph, sys.argv[2] + "_or")
    printGraphDatatable(compGraph, sys.argv[2] + "_co")
    printGraphDatatable(resGraph, sys.argv[2] + "_re")
else:
    printGraphDatatable(origGraph, "origGraph")
    printGraphDatatable(compGraph, "compGraph")
    printGraphDatatable(resGraph, "resGraph")


#print "Edgelist of original graph"
#print origGraph.get_edgelist()
#print "Edgelist of compacted graph"
#print compGraph.get_edgelist()
#print "Edgelist of resolved graph"
#print resGraph.get_edgelist()
#print "Problem Nodes in compacted graph:"
#print problemNodes
