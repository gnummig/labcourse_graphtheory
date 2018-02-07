#/bin/bash -x 

#find /homes/praktikumserver/praktikum/graph_prak_2018/ -name *.gtf 
echo header
while IFS= read -r GRAPH
do
	python parse3graphs.py $GRAPH 
done < $1

