#!/usr/bin/env python

import sys

f = open(sys.argv[1], 'r')
w = open(sys.argv[1]+".legend" , 'w')
while True:
	
	line = f.readline().strip() 
	if line == "@arcs":
		f.readline()
		break;

print "digraph g {"
print "graph [pad=\"0\", nodesep=\"0\", ranksep=\"0\"];"
print "rankdir=LR;"

i = 0
while True:
	raw = f.readline().strip()
 	if raw == '@attributes':
		break
	el = raw.split()
	if  el[5]=="4":
		print el[0],"->",el[1], ';'; 
	elif el[5]=="2":
		print el[0],"->",el[1], '[label="', el[6] ,'"];'; 
	else:	
		w.write(str(i)+"\t"+str(el[3][1:])+"\n")
		print el[0],"->",el[1], '[label="', "S"+str(i), el[6] ,'"];'; 
		i = i + 1
	

	
print "}"
