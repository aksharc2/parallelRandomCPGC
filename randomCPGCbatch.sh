#!/bin/bash
# Job name
#SBATCH --job-name parallel_random_CPGC
# Submit to the primary QoS
#SBATCH -q primary
# Request one node
#SBATCH -N 1
# Request entire node - the job allocation can not share nodes with other running jobs
#SBATCH --exclusive
# Total number of cores
#SBATCH -n 36
# Request memory
#SBATCH --mem=250G
# Mail when the job begins, ends, fails, requeues
#SBATCH --mail-type=ALL
# Where to send email alerts
#SBATCH --mail-user=gp6989@wayne.edu
# Create an output file that will be output_<jobid>.out
#SBATCH -o output_%j.out
# Create an error file that will be error_<jobid>.out
#SBATCH -e errors_%j.err
# Set maximum time limit
#SBATCH -t 72:00:00

cd /wsu/home/gp/gp69/gp6989/graphCompression/
#cd $TMPDIR

#cp -r /wsu/home/gp/gp69/gp6989/graphCompression/New_generated_data/ $TMPDIR
#cp /wsu/home/gp/gp69/gp6989/graphCompression/parallelCpaExe  $TMPDIR

echo "graphNodes,density,expNo,delta,compressionRatio,executionTime,readTime,writeTime,mergeTime,cores" > parallelRandomizedCPGCResults.csv
for exp in 1 2
do
    for node in 2048 4096 8192 16384 #32768
    do
        for density in 80 85 90 95 98 100 40 50 60 70
        do
            for run in 1 2 3 4 5 6 7 8 9 10
            do
                for delta in 0.5 0.6 0.7 0.8 0.9 1
                do
                    for proc in 4 8 12 16 20 24 28 32 36 # 40 44 48 52 56
                    do
                        fileName="datasets/bipartite_graph_${node}_${density}_${exp}.mtx"
                        mpirun -np $proc ./randomCPGC $fileName $node $delta $density  $exp >> parallelRandomizedCPGCResults.csv
			#sleep 1
                    done
                done
            done
        done
    done
done
