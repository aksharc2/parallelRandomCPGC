#!/bin/bash
# Job name
#SBATCH --job-name parallel_random_CPGC

# Submit to the primary QoS

#SBATCH -p RM

# Request one node

#SBATCH -N 1

# Total number of cores

#SBATCH -n 1

# Mail when the job begins, ends, fails, requeues

#SBATCH --mail-type=ALL

# Where to send email alerts

#SBATCH --mail-user=gw4590@wayne.edu

# Create an output file that will be output_<jobid>.out

#SBATCH -o output_%j.out

# Create an error file that will be error_<jobid>.out

#SBATCH -e errors_%j.err

# Set maximum time limit

#SBATCH -t 72:00:00

# Set minWall Time
#SBATCH --time-min=48:00:00

#cd /wsu/home/gw/gw45/gw4590/Graph_Compression/dataset1/
#cd /your current directory/Graph_Compression/datasets/
#echo "nodes,density,experimentNo"

set -x

for node in  32 64 128 #  256 512 1024 2048 4096 8192  32768 65536 131072 #262144
do
    for density in 80 85 90 95 98 #90 95  
    do
	for exp in 1 2 3 4 5 6 7 8 9 10
	do
	     python3 simpleGraphGenerator.py $node $density $exp
	done
    done
done
