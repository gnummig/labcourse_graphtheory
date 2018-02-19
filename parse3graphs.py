#!/usr/bin/env python2

import sys
from  networkx import *
import re
import os

#parse original graph
def getGraph():
    graph = MultiDiGraph()
    while True:
        line = f.readline().strip()
        if not line:
            break;
        if line == "@nodes":
            f.readline()
            break;
    while True:
        line = f.readline().strip()
        if not line : break;
        if line == "@arcs":
            f.readline()
            break;
        nodecount = int(line)
        graph.add_nodes_from(range(0, nodecount))

    while True:
        raw = f.readline().strip()
        if not line: break;
        if raw == '@attributes':
            break;
        el = raw.split()
        if  el[5]=="4":
            graph.add_edge( int( el[0] ) , int( el[1] ), type=4, label=el[2], Flow=0 )
        elif el[5]=="2":
            graph.add_edge( int( el[0] ) , int( el[1] ), type=2, label=el[2], Exon=el[3][1:], Flow=int( el[6]) )
            # todo does this work???
        else:
            #w.write(str(i)+"\t"+str(el[3][1:])+"\n")
            #print el[0],"->",el[1], '[label="', "S"+str(i), el[6] ,'"];';
            graph.add_edge( int( el[0] ) , int( el[1] ), type=1, label=el[2], binExon=el[3][1:], Flow=int( el[6]) )
    return graph

def detectProblemNodes(graph):
    inDegree = [b for (a,b) in list(graph.in_degree())]
    outDegree = [b for (a,b) in list(graph.out_degree())]
    problemNodes = [ 1*( a > 1 )*( b > 1 ) for a,b in zip( inDegree , outDegree ) ]
    for idx, isProblemNode in enumerate(problemNodes):
        graph.node[idx]['isProblemNode']=isProblemNode
    
def printGraphDatatable( graph , name ):
    inDegree = [b for (a,b) in list(graph.in_degree())]
    outDegree = [b for (a,b) in list(graph.out_degree())]
    problemNodes=get_node_attributes(graph,'isProblemNode')
    centrality=degree_centrality(graph).values()
    #GraphId, VertexCount, Vertex_ID, In_Degree, Out_Degree, ProblemNode?, Centrality
    for i in range( 0 , graph.number_of_nodes() ):
        print name +"\t"+ str( graph.number_of_nodes() ) + "\t" + str(i) + "\t" + str( inDegree[i] ) + "\t" + str( outDegree[i] ) + "\t" + str( problemNodes[i] ) + "\t" + str( centrality[i] )
        #print("")
    return

def createDotFile(graph, path):
    dotFile = open( path , 'w')
    dotFile.write("digraph g {\n")
    dotFile.write("graph [pad=\"0\", nodesep=\"0\", ranksep=\"0\"];\n")
    dotFile.write("rankdir=LR;\n")
    for (startnode, endnode , key) in list(graph.edges(keys=True)) :
        dotFile.write(str(startnode) + "->" + str(endnode) + '[label="' + str(graph[startnode][endnode][key]['Flow']) + '"];\n')
    isProblemNode = get_node_attributes(graph, 'isProblemNode').values();
    for idx  in range(0,graph.number_of_nodes()):
        if idx <= 1:
            dotFile.write('"' + str(idx) + '" [shape=circle, style=filled, fillcolor=blue]' )
        elif graph.nodes()[idx]['isProblemNode'] == 1:
            dotFile.write('"' + str(idx) + '" [shape=diamond, style=filled, fillcolor=red]' )
    dotFile.write("}")
    dotFile.close()

##########
#  main  #
##########


f = open( sys.argv[1] , 'r')
# remove first line, after that the exon list followx
f.readline().strip()
exonPos = []
while True:
    line = f.readline().strip()
    if not line:
        break;
    if "Bins" in line:
        break;
    exonPos.append( [ int( line.split(" ")[1].split("-")[0] ) , int( line.split(" ")[1].split("-")[1] ) ] )



origGraph = getGraph()
compGraph = getGraph()
resGraph = getGraph()

# get the transcripts
def getTranscripts( filename ):
    transcripts = open( os.path.dirname( sys.argv[1]) + filename , 'r')
    # structure of exonPos is a list with sublists [ [ exonpositions ] , FPKM ]
    ExonPos=[]
    transcriptcount=-1
    while True:
        line = transcripts.readline().strip()
        if not line:
            break;
        if "transcript\t" in line:
            line = transcripts.readline()
            transcriptcount = transcriptcount + 1
            ExonPos.append( [ [] , line.split(";")[2][7:-1] ] ) # FPKM value
        ExonPos[ transcriptcount ][0].append( int( line.split()[3] ) )
        ExonPos[ transcriptcount ][0].append( int( line.split()[4] ) )
        ExonPos[ transcriptcount ][0].sort()
    return ExonPos
transcriptExonPos = getTranscripts( "/transcripts.gtf" )
trueExonPos = getTranscripts( "/../truth.gtf" )
# todo neeed to flatten that list.
# translate the binary exon representation into actual exon posiions
spliceEdges = []
for (startnode, endnode , key) in list(resGraph.edges(keys=True)) :
    binex=resGraph[startnode][endnode][key]['binExon']
    #indeces of the exons on the path
    indices = [ b for a,b in zip( binex , range( 0 , len( binex ))) if int( a ) > 0 ]
    # [exonpos start,exonpos end]
    splicePosLong=[ b for a , b in zip( binex , exonPos ) if int( a ) > 0 ]
    # flatten list
    splicePosflat = [ int( item ) for sublist in splicePosLong for item in sublist ]
    # merge neighboring exons
    splicePosShort = [ a for a in splicePosflat if ( ( a + 1 ) not in splicePosflat ) & ( ( a - 1 ) not in splicePosflat ) ]
    # sort the list, since sometimes the true tanscripts are reversed
    splicePosShort.sort()
    # collect the edges of the resolved splicegraph
    spliceEdges.append( [ splicePosShort , [ startnode , endnode, resGraph[startnode][endnode][key]['label'], 0 ] ] )

def checktranscript( transcript ,transcriptPosition, graph , startnode , truePathVar ) :
    # if in the true trpts, take all paths, and mark used ones, if in the transpts, take only marked paths
    edgeset = [ edge for edge in  graph if edge[1][0] == startnode and ( edge[1][3] or truePathVar )  ]
    for edge in edgeset :
        thisedge = True # becomes false if one exon missmatches
        if len(transcript[0]) <= transcriptPosition + len(edge[0]) -2 :
            thisedge = False # if the edge is too long dont try to compare it, try the next edge
            continue;
        # if on the first edge, we need to check if the transcript really starts on the first edge
        if startnode==0:
            if transcript[0][0] > exonPos[0][1]:
                if not transcript[0][0]== edge[0][0]:
                    continue;
        for idx,splice in enumerate(edge[0][1:-1]): # check whether exons of the edge are in the transcript
            if not splice ==  transcript[0][transcriptPosition + idx +1]:
                thisedge = False # if any exon missmatches, try the next edge
                break; # dont try no more exons
        if not thisedge:
            continue; # with next edge
        if not edge[1][1] == 1: #' not at the end yet?
            # iteration 
            transcriptFound = checktranscript( transcript, transcriptPosition + len(edge[0][1:-1]) , graph , edge[1][1] , truePathVar )
            # if no path based on the current edge is found, continue and try another one
            if not transcriptFound:
                continue;# with next edge
            else:
                #append edge to list
                transcriptFound.append( edge[1][2] )
                edge[1][3]= edge[1][3] or truePathVar # if the true paths passes, this is a valid edge
                return transcriptFound
        # check if the transcript is covered til the end by this last exon
        elif not transcriptPosition + len(edge[0]) == len(transcript[0])  :
            continue;
        elif transcript[0][-1] < exonPos[-1][0]:
            if not transcript[0][-1] == edge[0][-1]:
                continue;
        else:
            # initialize with the edge that lead to endnode
            transcriptFound = [ edge[1][2] ]
            edge[1][3] = edge[1][3] or truePathVar
            return transcriptFound
    # there were no edges so this we wont get anywhere from here
    return False

def getTranscripPath( transcripts , truePathVar ):
    paths = []
    for transcript in transcripts:
        paths.append( checktranscript( transcript , 0 , spliceEdges , 0 , truePathVar ) )
    return paths

# takes a transcript path that is known to be chimaer, and a position from which ti start searching
def getChimaerNodes(transcriptPath , position):
    maxhops=0
    for truePath in truePaths:
        thispath=True
        if truePath:
            for idx,edge in enumerate(truePath[position:]):
                # search for the longest path you can go in one transcript from "position"
                if edge in transcriptPath:
                    maxhops = max(maxhops,idx+1)
                else:
                    break;
    # maxindex should now contain the maximum number of edges walked in one piece
    position = position + maxhops
    # if that does not take us to the end, iterate!
    if  position  != len(transcriptPath):
       chimaernodes = getChimaerNodes(transcriptPath, position)
       # find the end node of the last edge that matched, that is going to be the responsible one
       chimaernodes.append( [v for (u,v,c) in list(resGraph.edges(keys=True)) if resGraph[u][v][c]['label'] == transcriptPath [ maxhops ] ] )
    # if it is the end initialixe the list
    else:
        return []
    return chimaernodes


truePaths = getTranscripPath( trueExonPos , 1 )
transcriptPaths = getTranscripPath( transcriptExonPos, 0  )
print "true paths  through the graph"
print   truePaths
print "transcript paths that go only through edges of the graph that are used by truth"
print   transcriptPaths
chimaerNodesNested=[]
for idx,path in enumerate(transcriptPaths):
    if path and path in truePaths:
        print "correct transcript:"
        print path
    elif path:
        print "chimaer transcript:"
        print path
        chimaerNodesNested.append(getChimaerNodes(path , 0))
        print "the responsible nodes are:"
        # double flatten the list and remove duplicates by back and forth transforming to a set
        chimaerNodes = list(set([ node for sublist in chimaerNodesNested for subsublist in sublist for  node in subsublist] ) )
        print chimaerNodes
    else:
        print "incorrect transcript"
        print trueExonPos[idx]


detectProblemNodes(origGraph)
detectProblemNodes(compGraph)
detectProblemNodes(resGraph)

## should output a list of edges for each transcript and for the truth, how many transcripts have no path
    
graphIndex=""
for idx in range(0,len(sys.argv)-1):
    if sys.argv[idx] == "-dot":
        createDotFile(origGraph, sys.argv[idx+1] +  "_orig.dot")
        createDotFile(compGraph,  sys.argv[idx+1]  +  "_comp.dot")
        createDotFile(resGraph,  sys.argv[idx+1]  +  "_res.dot")
    if sys.argv[idx] == "-index":
        graphIndex = sys.argv[idx+1]
        
        
printGraphDatatable( origGraph , graphIndex + "orig" )
printGraphDatatable( compGraph , graphIndex + "comp" )
printGraphDatatable( resGraph , graphIndex + "res" )
        

#print "Edgelist of original graph"
#print origGraph.get_edgelist()
#print "Edgelist of compacted graph"
#print compGraph.get_edgelist()
#print "Edgelist of resolved graph"
#print resGraph.get_edgelist()
#print "Problem Nodes in compacted graph:"
#print problemNodes
# vim: noai:ts=4
