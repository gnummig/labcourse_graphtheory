#/bin/bash -x 

#find /homes/praktikumserver/praktikum/graph_prak_2018/ -name *.gtf 
echo "GraphID	 VertexCount	Vertex_ID	In_Degree	Out_Degree	ProblemNode	Centrality"
GraphID=0
while IFS= read -r GRAPH
do
	python parse3graphs.py $GRAPH -index $GraphID"_"
	let "GraphID++"
done < $1

