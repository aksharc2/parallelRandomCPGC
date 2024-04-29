# Paralle Randomized Clique Partition based Graph Compression (Par-CPGC)

## Introducation
The project aims to compress a given bipartiate graph by extracting cliqies. Doing so the compressed graph maintains the path infomation i.e., the connectivity of the vertices in the graph. This approach fall into lossless compression techniques. The compressed graph can be used as an input graph for other graph based algorithms such as matching, all pairs shortes path, etc. 
This git repository aim to provide the c scripts to obtain compression. It also has a python script for generating sample bipartite graphs.

Following are the programs in this repository:


1. randomCPGC.c (Our approach)
2. sequential_randomCPGC.c (sequentail version of our approach)
3. fm.c (baseline that is related to our approach presented in [1])


## Procedure to run the scripts
Instructions to run the script on a graph and obtain compression.

1. **Prerequisite packages and libraries required**
  - for c program
    - gcc version 9.x and above (the program is tested with gcc versions 11.4,  13.2)
    - OpenMPI version 4.0
  - for running python script
    - Numpy
    - Matplotlib.pyplot
    - Pandas
    - Statistics
    - Math
     

## References

[1] Tomás Feder and Rajeev Motwani. 1991. Clique partitions, graph compression and speeding-up algorithms. In Proceedings of the twenty-third annual ACM symposium on Theory of Computing (STOC '91). Association for Computing Machinery, New York, NY, USA, 123–133. [Link to paper](https://doi.org/10.1145/103418.103424).



# Parallel Randomized Algorithm for Graph Compression (Par-RCP)

## Abstract

This artifact aims to reproduce the experimental results presented in our paper, which proposes a new randomized parallel algorithm Par-RCP for graph compression based on partitioning the graph into bipartite cliques. We first design a randomized sequential algorithm Seq-RCP for graph compression based on clique partitioning and then use it as a base for the design of our parallel randomized algorithm Par-RCP. The algorithms extract δ-cliques C_q(U_q, W_q) from the bipartite graph G(U,W,E) and replace each of them with a tree formed by adding vertex z_q corresponding to the δ-clique C_q and connecting it to all the vertices in the two bipartitions of the clique. A δ-clique C_q in G is a complete bipartite subgraph with the left partition U_q of size ⌈n^(1-δ)⌉ and the right partition W_q of size k(n, m, δ) = ⌊δ log(n) / log(2n^2/m)⌋, where δ is a constant such that 0 ≤ δ ≤ 1, and m is the number of edges in G. We implement and perform an extensive performance evaluation on a system with a large number of cores. The results show that for large graphs with billions of edges, the parallel algorithm achieves a considerable speedup relative to its sequential counterpart of up to 29.25 when using 80 cores.

## Implementation and Experiments

### Artifact Dependencies and Requirements:

The artifact is tested on Linux systems with Ubuntu 20 or higher. For successfully running the programs, the machine would need a GCC compiler with version 9.4.0 or higher and an MPICC compiler with version 4.0 or higher.

#### Installing GCC compiler:

```bash
sudo apt update
sudo apt install build-essential
```

#### Installing OpenMPI:
```bash
sudo apt-get install openmpi-bin openmpi-doc libopenmpi-dev
```

### Generating Datasets:
A bipartite graph $G(U, W, E)$ can be generated using the Python script simpleGraphGenerator.py with three arguments in the following order:
- nodes: the number of vertices $n$, in left or right partition of the given graph.
- density: density $\rho$, i.e. the ratio of number of edges in given graph over the maximum number of possible edges ($n^2) in the given graph.
- experimentNo: it is an identifier of the generated graph with same nodes $n$ and density $\rho$.
```bash
python3 simpleGraphGenerator.py nodes density experimentNo
```

### Compiling and executing Seq-RCP, Par-RCP, and FM C programs Compiling Seq-RCP:
- Compiling Seq-RCP: ```bash gcc sequential_randomCPGC.c -o sequentialRCPGC -lm ```



























