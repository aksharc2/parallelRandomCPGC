# Paralle Randomized Clique Partition based Graph Compression (Par-CPGC)

## Introducation
The project aims to compress a given bipartiate graph by extracting cliqies. Doing so the compressed graph maintains the path infomation i.e., the connectivity of the vertices in the graph. This approach fall into lossless compression techniques. The compressed graph can be used as an input graph for other graph based algorithms such as matching, all pairs shortes path, etc. 
This git repository aim to provide the c scripts to obtain compression. It also has a python script for generating sample bipartite graphs.
Following are the programs in this repository:


1. randomCPGC.c (Our approach)
2. sequential_randomCPGC.c (sequentail version of our approach)
3. fm.c (baseline that is related to our approach presented in [1])




## References

[1] Tomás Feder and Rajeev Motwani. 1991. Clique partitions, graph compression and speeding-up algorithms. In Proceedings of the twenty-third annual ACM symposium on Theory of Computing (STOC '91). Association for Computing Machinery, New York, NY, USA, 123–133. [Link to paper](https://doi.org/10.1145/103418.103424).
