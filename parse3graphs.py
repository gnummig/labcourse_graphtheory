#!/usr/bin/env python2

import sys
from  igraph import *
import re

#parse original graph
def getGraph():
    graph = Graph( directed = True )
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
    graph.add_vertices( nodecount + 1 )
    edgecount = 0
    while True:
        raw = f.readline().strip()
        if not line: break;
        if raw == '@attributes':
            break;
        el = raw.split()
        graph.add_edges( [ (int( el[0] ) , int( el[1] ) ) ] )
        if  el[5]=="4":
            graph.es[ edgecount ]["type"] = 4
            graph.es[ edgecount ]["label"] = el[2]
            graph.es[ edgecount ]["Flow/Capacity"] = 0
        if el[5]=="2":
            graph.es[ edgecount ]["type"] = 2
            graph.es[ edgecount ]["label"] = el[2]
            graph.es[ edgecount ]["Exon"] = el[3][1:]
            graph.es[ edgecount ]["Flow/Capacity"] = int( el[6] )
            # todo does this work???
        else:
            #w.write(str(i)+"\t"+str(el[3][1:])+"\n")
            #print el[0],"->",el[1], '[label="', "S"+str(i), el[6] ,'"];';
            graph.es[ edgecount ]["type"] = 1
            graph.es[ edgecount ]["label"] = el[2]
            graph.es[ edgecount ]["Flow/Capacity"] = int(el[6])
            graph.es[ edgecount ]["binExon"] = el[3][1:]
        edgecount = edgecount + 1
    return graph

def printGraphDatatable( graph , name ):
    inDegree =  graph.degree( type = "in" )
    outDegree = graph.degree( type = "out" )
    problemNodes = [ 1*( a > 1 )*( b > 1 ) for a,b in zip( inDegree , outDegree ) ]
    centrality = graph.evcent( directed = False , weights = "Flow/Capacity" )

    #GraphId, VertexCount, Vertex_ID, In_Degree, Out_Degree, ProblemNode?, Centrality
    for i in range( 0 , graph.vcount() ):
        foo = ""
#print name +"\t"+ str( graph.vcount() ) + "\t" + str(i) + "\t" + str( inDegree[i] ) + "\t" + str( outDegree[i] ) + "\t" + str( problemNodes[i] ) + "\t" + str( centrality[i] )
        #print("")
    return

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
    return ExonPos
transcriptExonPos = getTranscripts( "/transcripts.gtf" )
trueExonPos = getTranscripts( "/../truth.gtf" )
# todo neeed to flatten that list.
# translate the binary exon representation into actual exon posiions
spliceEdges = []
for idx , binex in enumerate( resGraph.es[ "binExon" ] ) :
    #indeces of the exons on the path
    indices = [ b for a,b in zip( binex , range( 0 , len( binex ))) if int( a ) > 0 ]
    # [exonpos start,exonpos end]
    splicePosLong=[ b for a , b in zip( binex , exonPos ) if int( a ) > 0 ]
    # flatten list
    splicePosflat = [ int( item ) for sublist in splicePosLong for item in sublist ]
    # merge neighboring exons
    splicePosShort = [ a for a in splicePosflat if ( ( a + 1 ) not in splicePosflat ) & ( ( a - 1 ) not in splicePosflat ) ]
    # collect the edges of the resolved splicegraph
    spliceEdges.append( [ splicePosShort , [ resGraph.es[idx].source , resGraph.es[idx].target,resGraph.es[idx]["label"] , 0 ] ] )

def checktranscript( transcript ,transcriptPosition, graph , startnode , truePathVar ) :
    # if in the true trpts, take all paths, and mark used ones, if in the transpts, take only marked paths
    edgeset = [ edge for edge in  graph if edge[1][0] == startnode and ( edge[1][3] or truePathVar )  ]
    for edge in edgeset :
        thisedge = True # becomes false if one exon missmatches
        if len(transcript[0]) <= transcriptPosition + len(edge[0]) -2 :
            thisedge = False # if the edge is too long dont try to compare it, try the next edge
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

truePaths = getTranscripPath( trueExonPos , 1 )
transcriptPaths = getTranscripPath( transcriptExonPos, 0  )
print "transcript paths through the graph"
print   transcriptPaths
print "true paths through the graph"
print   truePaths
for idx,path in enumerate(transcriptPaths):
    if path and path in truePaths:
        print "correct transcript:"
        print path
    elif path:
        print "chimaer transcript:"
        print path
    else:
        print "incorrect transcript"
        print trueExonPos[idx]




## should output a list of edges for each transcript and for the truth, how many transcripts have no path

if  len( sys.argv ) > 2 :
    printGraphDatatable( origGraph , sys.argv[2] + "_or" )
    printGraphDatatable( compGraph , sys.argv[2] + "_co" )
    printGraphDatatable( resGraph , sys.argv[2] + "_re" )
else:
    printGraphDatatable( origGraph , "origGraph" )
    printGraphDatatable( compGraph , "compGraph" )
    printGraphDatatable( resGraph , "resGraph" )


#print "Edgelist of original graph"
#print origGraph.get_edgelist()
#print "Edgelist of compacted graph"
#print compGraph.get_edgelist()
#print "Edgelist of resolved graph"
#print resGraph.get_edgelist()
#print "Problem Nodes in compacted graph:"
#print problemNodes
# vim: noai:ts=4
