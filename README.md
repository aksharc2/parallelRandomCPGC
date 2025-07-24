# Parallel Randomized Clique Partitioning-Based Algorithm for Graph Restructuring (Par-CPGC)

## Abstract

This artifact aims to reproduce the experimental results presented in our paper, which proposes a new randomized parallel algorithm Par-RCP for graph restructuring based on partitioning the graph into bipartite cliques. We first design a randomized sequential algorithm Seq-RCP for graph restructuring based on clique partitioning and then use it as a base for the design of our parallel randomized algorithm Par-RCP. The algorithms extract δ-cliques C_q(U_q, W_q) from the bipartite graph G(U,W,E) and replace each of them with a tree formed by adding vertex z_q corresponding to the δ-clique C_q and connecting it to all the vertices in the two bipartitions of the clique. A δ-clique C_q in G is a complete bipartite subgraph with the left partition U_q of size ⌈n^(1-δ)⌉ and the right partition W_q of size k(n, m, δ) = ⌊δ log(n) / log(2n^2/m)⌋, where δ is a constant such that 0 ≤ δ ≤ 1, and m is the number of edges in G. Extensive performance evaluation on  large-scale graphs with billions of edges shows that \textsf{Par-RCP} reduces the edge count by up to 6.69$\times$ and achieves a speedup of up to 29.25. Applying the restructured graph further improves algorithmic performance, with Dinitz's algorithm achieving up to 2.06 speedup and All-Pairs Shortest Paths up to 213.34 relative to the case of using the original input graph as input. These results highlight the scalability and practical utility of our approach for large-scale graph analytics.


## Introduction
The project aims to restructure a given bipartiate graph by extracting bipartite cliqies and replacing them with tripartite structure. Doing so the restructured graph maintains the path infomation i.e., the connectivity of the vertices in the graph. This approach fall into lossless restructuring techniques. The restructured graph can be used as an input graph for other graph based algorithms such as matching, all pairs shortes path, etc. 
This git repository aim to provide the c scripts to obtain the restrectured graph. It also has a python script for generating sample bipartite graphs.

Following are the programs in this repository:
1. randomCPGC.c (Our approach)
2. sequential_randomCPGC.c (sequentail version of our approach)
3. fm.c (baseline that is related to our approach presented in [1])


## Procedure to run the scripts
Instructions to run the script on a given original graph and obtain the restrectured graph.

1. **Prerequisite packages and libraries required**
  - Tested with Ubuntu 20.04 and above versions
  - for c program
    - gcc version 9.x and above
    - OpenMPI version 4.0
  - for running python script (used for generating plots for the paper)
    - Numpy
    - Matplotlib.pyplot
    - Pandas
    - Statistics
    - Math
     

## Implementation and Experiments

### Artifact Dependencies and Requirements:

The artifact is tested on Linux systems with Ubuntu 20 or higher. For successfully running the programs, the machine would need a GCC compiler with version 9.4.0 or higher and an MPICC compiler with version 4.0 or higher.

#### Installing GCC compiler:
- ```bash
  sudo apt update
  sudo apt install build-essential
  ```
- For detailed installation procedure see [How to Install GCC Compiler on Ubuntu?](https://phoenixnap.com/kb/install-gcc-ubuntu)

#### Installing OpenMPI:
- ```bash
  sudo apt-get install openmpi-bin openmpi-doc libopenmpi-dev
  ```
- For detailed installation procedure see [Installing OpenMPI](https://webpages.charlotte.edu/abw/coit-grid01.uncc.edu/ParallelProgSoftware/Software/OpenMPIInstall.pdf)

### Generating Datasets:
A bipartite graph $G(U, W, E)$ can be generated using the Python script simpleGraphGenerator.py with three arguments in the following order:
- nodes: the number of vertices $n$, in left or right partition of the given graph.
- density: density $\rho$, i.e. the ratio of number of edges in given graph over the maximum number of possible edges ($n^2) in the given graph.
- experimentNo: it is an identifier of the generated graph with same nodes $n$ and density $\rho$.
```bash
python3 simpleGraphGenerator.py nodes density experimentNo
```

### Compiling and executing Seq-RCP, Par-RCP, and FM C programs Compiling Seq-RCP:
- Compiling Seq-RCP: ```gcc sequential_randomCPGC.c -o sequentialRCPGC -lm ```
- Executing Seq-RCP: ```./sequentialRCPGC fileName node delta density exp ```
  - fileName: is the name of the input file along with the absolute or relative path.
  - nodes: the number of vertices $n$, in left or right partition of the given graph.
  - delta: it is constant factor $\delta$, i.e. $0 < \delta < 1$.
  - density: density $\rho$, i.e. the ratio of number of edges in given graph over the maximum number of possible edges ($n^2$) in the given graph.
  - exp: the instance of the generated graph with same nodes $n$ and density $\rho$.

- Compiling Par-RCP: ```mpicc randomCPGC.c -o randomCPGC -lm ```
- Executing Par-RCP: ```mpirun -np nproc ./randomCPGC fileName node delta density exp ```
  - nproc: number of processors to use for running the program.
  - fileName: is the name of the input file along with the absolute or relative path.
  - nodes: the number of vertices $n$, in left or right partition of the given graph.
  - density: density $\rho$, i.e. the ratio of number of edges in given graph over the maximum number of possible edges ($n^2$) in the given graph.
  - exp: the instance of the generated graph with same nodes $n$ and density $\rho$.

- Compiling FM: ```gcc fm.c -o fm -lm```
- Executing FM: ```./fm fileName node delta density exp```
  - fileName: is the name of the input file along with the absolute or relative path.
  - nodes: the number of vertices $n$, in left or right partition of the given graph.
  - density: density $\rho$, i.e. the ratio of number of edges in given graph over the maximum number of possible edges ($n^2$) in the given graph.
  - exp: the instance of the generated graph with same nodes $n$ and density $\rho$.

### Compiling and executing Dinitz and All-Pairs Shortest Path (APSP) algorithms:
- Compiling Dinitz for bipartite graph: ```gcc dinics_bi.c -o dinitz_bi -lm ```
- Compiling Dinitz for tripartite graph: ```gcc dinics_tri.c -o dinitz_tri -lm ```
- Executing Dinitz for bipartite graph: ```./dinitz_bi node density exp delta ```
- Executing Dinitz for tripartite graph: ```./dinitz_tri node density exp delta ```  ./c_apsp "$file" -d 1
  - nodes: the number of vertices $n$, in left or right partition of the given graph.
  - delta: it is constant factor $\delta$, i.e. $0 < \delta < 1$.
  - density: density $\rho$, i.e. the ratio of number of edges in given graph over the maximum number of possible edges ($n^2$) in the given graph.
  - exp: the instance of the generated graph with same nodes $n$ and density $\rho$.
- Compiling APSP: ```gcc apsp_bsf_reconstructed_graph.c -lm -o apsp ```
- Executing APSP: ```./apsp fileName -d 1 -c 1 ```
  - fileName: is the name of the input file along with the absolute or relative path.
  - directed (-d) option: use '1' for using directed graph as input and '0' for non-directed graph.
  - restructured (-c) option: use '1' for using restructed graph as input and '0' for given original graph.
 
    
### Running multiple experiments
Following Bash scripts allows to run multiple experiments with Seq-RCP, Par-RCP, and FM and saves the results in csv files.

- Running multiple experiments with Seq-RCP and Par-RCP:
  - use the bash script "randomCPGCbatch.sh"
  - update the for loops for experimentno, nodes, density, delta and number of processors to use (only for running Par-RCP) as shown in the image below.
  - run the script with ```bash randomCPGCbatch.sh``` command
  - when this program termintaes it creates following files:
    - **parCPGC_results.csv** which has the execution time and the restructured graph obtained for particular instance of the experiments from the parallel program.
    - restructured graphs from the parallel programs are stored in the restrecturedGraph directory in this git repository directory and are named with respective instance of the given graph, for example ***restructured_graph_32_80_1_0.9.mtx***.
    - **seqCPGC_results.csv** which has the execution time and the restructured graph obtained for particular instance of the experiments from the sequential program.
    - restructured graphs from the sequentaial programs are stored in the git repository directory and are named with respective instance of the given graph, for example ***S_restrectured_graph_32_80_1.mtx***.
   
- Running multiple experiments with FM:
  - use the bash script "fmbatchScript.sh"
  - update the for loops for experimentno, nodes, density, delta and number of processors to use (only for running Par-RCP) as shown in the image below.
  - run the script with ```bash fmbatchScript.sh``` command
  - when this program termintaes it creates following files:
    - **fm_results.csv** which has the execution time and the restructured graph obtained for particular instance of the experiments from the parallel program.
    - restructured graphs from the parallel programs are stored in the **fm_restructured_graphs** directory in this git repository directory and are named with respective instance of the given graph, for example ***tripartite_graph_32_80_1_90.mtx***.

### Results
The obtained edge reduction ratio and the running time of the C programs for Seq-RCP, Par-RCP, and FM are stored in ".csv" files named "seqCPGC_results.csv", "parCPGC_results.csv", and "fm_results.csv", respectively.


## Generating plots with the results
We use Jupyter Notebook to generate plots for the paper. Following instrustions provide Jupyter Notebook installation procedure and how to generate the plots:
- installing Jupyter Notebook 
  - ```bash
    sudo apt-get update
    sudo apt-get install python3
    sudo apt-get install python3-pip
    pip3 install notebook
    ```
  - If you see a warning related to path as shown in Image: 3, then add the path where the package is installed to the PATH using the following commands
    - ```bash export PATH="${PATH}:/home/cpgc/.local/bin"  ```, replace "/home/cpgc/.local/bin" with the path where Jupyter Notebook is installed.
    - ```bash source ~/.bashrc```
  - For detailed installation procedure see [Jupyter Notebook Installation Procedure](https://saturncloud.io/blog/how-to-install-jupyter-notebook-in-ubuntu/)
- The jupyter notebook files have command for installing the required packages for generating the plots. Incase the package installation is unsucessful then follow the following procedure to install respective package:
  - **Pandas installation:**
    - ```bash pip install pandas```
    - For detailed installation procedure see [Link](https://pandas.pydata.org/docs/getting_started/install.html)
  - **Numpy installation:**
    - ```bash pip install numpy```
    - For detailed installation procedure see [Link](https://numpy.org/install/)
  - **Matplotlib installation:**
    - ```bash pip install matplotlib```
    - For detailed installation procedure see [Link](https://pypi.org/project/matplotlib/)
  - **Statistics installation:**
    - ```bash pip install statistics```
    - For detailed installation procedure see [Link](https://pypi.org/project/statistics/)
  - **Math installation:**
    - ```bash pip install python-math```
    - For detailed installation procedure see [Link](https://pypi.org/project/python-math/)
- To generate results
  - open Jupyter Notebook by running this command in the terminal ```bash jupyter notebook``` which will open a web browser for Jupyter Notebook with /home directory.
  - On the opened web page with /home directory go to "parallelRandomCPGC-main" (or the directory when this git repo is saved/loaded) directory. If results are generated and stored in the respective files and location mentioned above, then the following _.ipynb_ files execution should present the result plots.
  - To generate result plots for Par-RCP open _Par_RCPGC_plots.ipynb_ file and run all cells.
  - To generate results comparing the baseline FM open _comparison_with_fm.ipynb_ file and run all cells.
  - The _get_k_hat.ipynb_ plots the $\hat{k}$ comparision of our approach with FM for given nodes $n$, edges $m$, density $\rho$, and $\delta$.
    - The function **get_k()** provides a starting value of $\hat{k}$.



## References

[1] Tomás Feder and Rajeev Motwani. 1991. Clique partitions, graph compression and speeding-up algorithms. In Proceedings of the twenty-third annual ACM symposium on Theory of Computing (STOC '91). Association for Computing Machinery, New York, NY, USA, 123–133. [Link to paper](https://doi.org/10.1145/103418.103424).

















