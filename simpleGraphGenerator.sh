#cd /wsu/home/gw/gw45/gw4590/Graph_Compression/dataset1/
#cd /your current directory/Graph_Compression/datasets/
#echo "nodes,density,experimentNo"

for node in  16384 32768 #32 64 128 256 512 1024 2048 4096 8192
do
    for density in  80 85 90 95 98
    do
	for exp in  1 2 3 4 5 6 7 8 9 10
	do
	     python3 simpleGraphGenerator.py $node $density $exp
	done
    done
done
