#/bin/bash -x 

#find /homes/praktikumserver/praktikum/graph_prak_2018/ -name *.gtf 
echo "GraphID	GraphKind	VertexCount	VertexID	ProblemNode	ChimearNode	In_Degree	Out_Degree	In_Flow	In_Flow_Std	Out_Flow	Out_Flow_Std	Centrality"
GraphID=0
while IFS= read -r GRAPH
do
	python parse3graphs.py $GRAPH -index $GraphID
	let "GraphID++"
done < $1

