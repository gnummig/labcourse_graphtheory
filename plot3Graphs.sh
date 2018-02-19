#/bin/bash -x 

python parse3graphs.py $1 -dot $1
dot $1"_orig.dot" -Tpdf > $1"_orig.pdf" 
dot $1_"comp.dot" -Tpdf > $1_"comp.pdf"
dot $1_"res.dot" -Tpdf > $1_"res.pdf"
rm $1"_orig.dot" $1_"comp.dot" $1_"res.dot"
convert -background none -gravity center  -density 300  $1"_orig.pdf" $1_"comp.pdf" $1_"res.pdf" -append -trim $1".png"