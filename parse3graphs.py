#!/usr/bin/env python2

import sys
from  networkx import *
import numpy
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
        nodecount = int(line) + 1
        graph.add_nodes_from(range(0, nodecount), isChimearNode=0, isProblemNode=0)

    while True:
        raw = f.readline().strip()
        if not line: break;
        if raw == '@attributes':
            break;
        el = raw.split()
        if  el[5]=="4":
            Exon = bin( 0 )[ 2: ].zfill( len( exonPos ) )
            graph.add_edge( int( el[0] ) , int( el[1] ), type=4, label=el[2], binExon = Exon, Flow=0 )
        elif el[5]=="2":
            Exon = bin( int( el[4][:-1] ) )[2:].zfill( len( exonPos ) )
            graph.add_edge( int( el[0] ) , int( el[1] ), type=2, label=el[2], binExon = Exon, Flow=int( el[6]) )
            # todo does this work???
        else:
            #w.write(str(i)+"\t"+str(el[3][1:])+"\n")
            #print el[0],"->",el[1], '[label="', "S"+str(i), el[6] ,'"];';
            graph.add_edge( int( el[0] ) , int( el[1] ), type=1, label=el[2], binExon=el[3][1:], Flow=int( el[6]) )
    #delete nodes without edges:
    graph.remove_nodes_from([u for u,v in graph.degree() if v==0 ])
    return graph

def computeGraphAttributes(graph):
    for v in graph.nodes():
        graph.nodes()[v]['isProblemNode'] = 1*( graph.in_degree()[v] > 1 )*( graph.out_degree()[v]  > 1 )
        outFlows=[f for u1, u2, f in graph.out_edges(v, data='Flow')]
        inFlows=[f for u1, u2, f in graph.in_edges(v, data='Flow')]
        graph.nodes()[v]['inFlow'] = sum(inFlows)
        graph.nodes()[v]['outFlow'] = sum(outFlows)
        if inFlows:
            graph.nodes()[v]['inFlowStd'] = numpy.std(inFlows)
        # todo empty list is also true ?!?
        else:
            graph.nodes()[v]['inFlowStd'] = "NA"
        if outFlows:
            graph.nodes()[v]['outFlowStd'] = numpy.std(outFlows)
        else:
            graph.nodes()[v]['outFlowStd'] = "NA"


def printGraphDatatable( graph , name, graphKind):
    inDegree = [b for (a,b) in list(graph.in_degree())]
    outDegree = [b for (a,b) in list(graph.out_degree())]
    problemNodes=get_node_attributes(graph,'isProblemNode')
    chimearNodes=get_node_attributes(graph,'isChimearNode')
    centrality=degree_centrality(graph)
    inFlow = get_node_attributes(graph,'inFlow')
    inFlowStd = get_node_attributes(graph,'inFlowStd')
    outFlow = get_node_attributes(graph,'outFlow')
    outFlowStd = get_node_attributes(graph,'outFlowStd')
    #GraphId, GraphKind , VertexCount, VertexID, ProblemNode?, ChimearNode? , In_Degree, Out_Degree, In_Flow, In_Flow_Std, Out_Flow, Out_Flow_Std, Centrality
    for v in graph.nodes():
        print name + "\t" + graphKind + "\t" + str( graph.number_of_nodes() ) + "\t" + str(v) +"\t" + str( problemNodes[v] ) + "\t" + str( chimearNodes[v] ) + "\t" + str( graph.in_degree[v] ) + "\t" + str( graph.out_degree[v] ) + "\t"  + str(inFlow[v]) + "\t"  + str(inFlowStd[v]) + "\t"  + str(outFlow[v]) + "\t"  + str(outFlowStd[v]) + "\t" + str( centrality[v] )
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
        elif graph.nodes()[idx]['isChimearNode'] > 0:
            dotFile.write('"' + str(idx) + '" [shape=diamond, style=filled, fillcolor=red]' )
        elif graph.nodes()[idx]['isProblemNode'] == 1:
            dotFile.write('"' + str(idx) + '" [shape=diamond, style=filled, fillcolor=orange]' )
    dotFile.write("}")
    dotFile.close()

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

## should output a list of edges for each transcript and for the truth, how many transcripts have no path
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

#basically just a function to parse transcripts to the  checktranscript function
def getTranscripPath( transcripts , truePathVar ):
    paths = []
    for transcript in transcripts:
        paths.append( checktranscript( transcript , 0 , spliceEdges , 0 , truePathVar ) )
    return paths

# takes a transcript path that is known to be chimaer, and a position from which to start searching
def getChimaerNodes(chimaerPath ):
#tuple [truepath just walked, chimaer nodes] 
    chimaerwalk = [ [ -1 , [] ] ]
    for edge in chimaerPath :
        justwalked = []
        for true_path_num , truePath in enumerate(truePaths):
            # for every truepath that we can go calclute the best paths 
            if edge in truePath:
                possibletrue = []
                # construct list of transitions from  oldpaths to this truePath
                for  [oldpath,nodeslist] in chimaerwalk :
                    if ( oldpath == true_path_num ) or ( oldpath == -1 ) :
                        possibletrue.append( [ true_path_num , nodeslist ] )
                    else:
                        # get the node that the chimaer edge came from (in reverse)
                        [new_cnode] =  [ v for (u,v,c) in list(resGraph.edges(keys=True)) if resGraph[u][v][c]['label'] == edge ]
                        newnodeslist =  list(nodeslist)
                        newnodeslist.append(new_cnode)
                        possibletrue.append( [true_path_num , newnodeslist] )
                bestscore = min( [ len(bb) for [aa,bb] in possibletrue ] )
                # from all combinations that end on this truepath take those with minimal cost
                for [d,e] in possibletrue:
                    if len(e) == bestscore:
                        justwalked.append( [d,e] )
        chimaerwalk = list(justwalked)
    final = []
    bestscore = min( [ len(b) for [a,b] in chimaerwalk ] )
    for [ a , b ] in chimaerwalk :
        if len( b ) == bestscore :
            final.append( b )
    return final
##########
#  main  #
##########
# chr2_0.graph

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

# get the transcripts that are there and that have been found
transcriptExonPos = getTranscripts( "/transcripts.gtf" )
trueExonPos = getTranscripts( "/../truth.gtf" )

# translate the binary exon representation of the resolved graph into actual exon posiions
spliceEdges = []
for (startnode , endnode , key ) in list( resGraph.edges( keys = True ) ) :
    binex=resGraph[ startnode ] [ endnode ][ key ] [ 'binExon' ]
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
# function to get the initial nodes that correspond to problemnodes in the final graph

def matchResolvedEdge( startnode , res_bin_edge , merged_bin_init_edge ) :
    for (startnode , endnode , key ) in list( origGraph.edges( keys = True ) ) :
        # add new exons to the edge we construct
        this_edge = [ a | b for a,b in zip( merged_bin_init_edge , origGraph[ startnode ][ endnode ][ key ][ 'binExon' ] ) ]
        # check if the strings match up to the last 1 in this_edge
        if [ a ^ b for a,b in zip( res_edge , this_edge ) ].index( 1 ) > len( this_edge ) - this_edge[::-1].index(1):
            # go for the next edge, as this one was fine.
            if this_edge == res_edge :  # wonder if that is correct, ie whether both are just arrays, but they should
                print, "that wonder worked, the edges matched"
                # initialize for succesful backward recursion
                origGraph.edges()[ ( startnode , endnode , key ) ] [ "inResolved" ] = 1
                return True
            #if we are not at the end of resedge, keep going
            result = matchResolvedEdge( endnote , res_edge , this_edge )
            if result:
                origGraph.edges()[ ( startnode , endnode , key ) ] [ "inResolved" ] = 1
                return True
    return False

#initialize new property " is in resolved" as false
for (startnode , endnode , key ) in list( origGraph.edges( keys = True ) ) :
    origGraph.edges()[ ( startnode , endnode , key ) ] [ "inResolved" ] = 0
# initialize a binary 0 string of the length of the others
# find all edges in initial graph that are compacte onto resolved ones, ie those who are not deleted, during resolvation
for (res_startnode , res_endnode , res_key ) in list( resGraph.edges( keys = True ) ) :
    # we surely find the startnode in res in the initial graph, and thats were our path starts
    matchResolvedEdges(res_startnode, (res_startnode , res_endnode , res_key),

# mark initial nodes that are compacted NODES that are compacted onto resolved with the node they become
for (startnode , endnode , key ) in [ (startnode , endnode , key ) for ( startnode , endnode , key ) in list( origGraph.edges( keys = True ) ) if  origGraph.edges()[ ( startnode , endnode , key ) ] [ "inResolved" ] == 1 ]:
    



truePaths = getTranscripPath( trueExonPos , 1 )
truePaths = [path for path in truePaths if path]
transcriptPaths = getTranscripPath( transcriptExonPos, 0  )
#print "true paths  through the graph"
#print   truePaths
#print "transcript paths that go only through edges of the graph that are used by truth"
#print   transcriptPaths
chimaerNodesNested=[]
for idx,path in enumerate(transcriptPaths):
    if path and path in truePaths:
        #print "correct transcript:"
        #print path
        continue;
    elif path:
        #print "chimaer transcript:"
        #print path
        chimaerNodesNested.append(getChimaerNodes(path))
        #print chimaerNodesNested
        #print "the responsible nodes are:"
        # double flatten the list and remove duplicates by back and forth transforming to a set
        chimaerNodes = list(set([ node for sublist in chimaerNodesNested for subsublist in sublist for  node in subsublist] ) )
        for v in chimaerNodes:
            resGraph.nodes()[v]['isChimearNode']=resGraph.nodes()[v]['isChimearNode']+1
        #print chimaerNodes
    #else:
        #print "incorrect transcript"
        #print trueExonPos[idx]


computeGraphAttributes(origGraph)
computeGraphAttributes(compGraph)
computeGraphAttributes(resGraph)


## should output a list of edges for each transcript and for the truth, how many transcripts have no path
graphIndex=""
for idx in range(0,len(sys.argv)-1):
    if sys.argv[idx] == "-dot":
        createDotFile(origGraph, sys.argv[idx+1] +  "_orig.dot")
        createDotFile(compGraph,  sys.argv[idx+1]  +  "_comp.dot")
        createDotFile(resGraph,  sys.argv[idx+1]  +  "_res.dot")
    if sys.argv[idx] == "-index":
        graphIndex = sys.argv[idx+1]


printGraphDatatable( origGraph , graphIndex, "orig" )
printGraphDatatable( compGraph , graphIndex, "comp" )
printGraphDatatable( resGraph , graphIndex, "res" )


#print "Edgelist of original graph"
#print origGraph.get_edgelist()
#print "Edgelist of compacted graph"
#print compGraph.get_edgelist()
#print "Edgelist of resolved graph"
#print resGraph.get_edgelist()
#print "Problem Nodes in compacted graph:"
#print problemNodes
# vim: noai:ts=4
