/*On windows follow this instruction to get a running time:
  1.clock_t start = clock();
  2.clock_t stop = clock();
  3.elapsed = ((double)(stop - start)) / CLOCKS_PER_SEC * 1000.0;
 */

 /*On Grid follow this instruction to get a running time:
   0. struct timespec begin, end;
   1.clock_gettime(CLOCK_REALTIME, &begin);
   2.clock_gettime(CLOCK_REALTIME, &end);         
   3.long seconds = end.tv_sec - begin.tv_sec;
   4.long nanoseconds = end.tv_nsec - begin.tv_nsec;
   5.elapsed = (seconds + nanoseconds * 1e-9) * 1000;
  */

#pragma warning(disable:4996)
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <stdbool.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <string.h>
#include <dirent.h>
#include <time.h>


#define MAXCHAR = 256

void createFolder(const char* folderName) {
    // Check if the folder exists, create it if not
    struct stat st;
    if (stat(folderName, &st) == -1) {
		mkdir(folderName, 0777);
    }
}


int get_k_hat(int graph_nodes, int m_hat, float delta) {
	int k_hat;
    float de = (2 * pow((double)graph_nodes, 2)) / m_hat;
    float nu = delta * log2((double)graph_nodes);
    k_hat = floor((double)nu / (log2((double)de)));
	return k_hat;
}


bool** adj_matrix;
int* leftClique;
int* rightClique;
int k_hat;
int delta;
int m_hat = 0;
int initialEdges;
int edgesRemoved = 0;
int cliqueEdges = 0;
int q; // Clique index
int* S;

int graphNodes;
struct timespec begin, end;
char f_name[100]; // name of the adjacency matrix file


void load_adj_matrix(){
    FILE* filePointer = fopen(f_name, "r");
    char line[MAXCHAR];
	int lp, rp, mp, edges, u, w;
    if (!filePointer)
        printf("Can't open the file\n");
    else {
        while (fgets(line, sizeof(line), filePointer)) {
			if (line[0] != '%'){
				break;
			}
        }
		sscanf(line, "%d %d %d", &lp, &rp, &edges);
		if (rp == 0)
			sscanf(line, "%d %d %d %d", &lp, &mp, &rp, &edges);
		for(int i = 0; i < edges; i++){
			fgets(line, sizeof(line), filePointer);
			sscanf(line, "%d %d", &u, &w);
			adj_matrix[u-1][w-1] = 1;
		}
		initialEdges = m_hat = edges;
        fclose(filePointer);
    }
}

void get_k_hat() {
    float de = (2 * pow((double)graphNodes, 2)) / m_hat;
    float nu = delta * log2((double)graphNodes);
    k_hat = floor((double)nu / (log2((double)de)));
}


bool** getAllocate(int n) {
    int i;
    bool** arr = (bool**)malloc((n) * sizeof(bool*));
    for (i = 0; i < n; i++) {
        arr[i] = (bool*)malloc((n) * sizeof(bool));
    }
    return arr;
}


void getDeAllocate(int n, bool** arr) {
    int i;
    for (i = 0; i < n; i++) {
        free(arr[i]);
    }
    free(arr);
}

int findCommonNeighbours() {
    int i, j;
    printf("\nLeft Partition: ");
    int u = 0;
    for (i = 0; i < graphNodes; i++) {
        int flag = 0;
        for (j = 0; j < k_hat && !flag; j++) {
            if (adj_matrix[i][rightClique[j]] == 0) {
                flag = 1;
            }
        }
        if (!flag) {
            leftClique[u] = i;
            for (j = 0; j < k_hat; j++) {
                adj_matrix[i][rightClique[j]] = 0;
                m_hat -= 1; // need to remove this for parallel implementation
                edgesRemoved += 1; // need to remove this for parallel implementation
            }
            printf("%d ", leftClique[u]);
            u++;
        }
    }
	cliqueEdges += u * k_hat;
    //printf("\n Size of Left Partition of %d Clique: %d", d, u);
    //printf("\nEdges in Clique: %d", cliqueEdges);
	return u;
}

float compressionRatio() {
    return ((float)initialEdges / (m_hat + cliqueEdges));
}

void saveCliqueEdges(int commonNeighbours, int q){
	for(int u = 0; u < commonNeighbours; u++){
		fprintf("%d %d\n", leftClique[u], q);
	}
	for (int w = 0; w < k_hat; w++){
		fprintf("%d %d\n", q, rightClique[w]);
	}
}


void RandomizedAlgorithm() {
    int nodeIndex;
	int q = graphNodes;
	get_k_hat();
    while (k_hat > 1) {
		for (int i = 0; i < graphNodes; i++){
			S[i] = i;
		}
		int setSize = graphNodes;
		while (setSize > 0){
			printf("Clique %d \nright partition: ", q);
			leftClique = (int*)malloc(graphNodes * sizeof(int));
			rightClique = (int*)malloc(k_hat * sizeof(int));
			if (setSize > k_hat){
				for(int i = 0; i < k_hat; i++){
					nodeIndex = rand() % setSize
					rightClique[i] = S[nodeIndex];
					setSize--;
					for (int j = nodeIndex; j < setSize; j++){
						S[j] = S[j+1];
					}
				}
			}
			else{
				for(int i = 0; i < setSize; i++){
					rightClique[i] = S[i];
					setSize--;
				}
			}
			int commonNeighbours = findCommonNeighbours();
			get_k_hat();
			saveCliqueEdges(commonNeighbours, q);
			free(leftClique);
			free(rightClique);
			q++;
		}
    }
}


int main() {
    int nodes = 32; // atoi(argv[1]);  // int argc, char* argv[]
    int density = 80; // atoi(argv[2]);
    int exp = 1; // atoi(argv[3]);
    srand(time(NULL));
    graphNodes = nodes;
    adj_matrix = getAllocate(graphNodes); 
    S = (int*)malloc((graphNodes) * sizeof(int));
    FILE* tempFile = fopen(filePath, "a");

    //sprintf(f_name, "Bipartite%dX%d.csv", graphNodes, graphNodes);
    sprintf(f_name, "New_generated_data/Bipartite_%dX%d/%d/Bipartite_%dX%d_%d_%d.csv", nodes, nodes, density, nodes, nodes, density, exp);

    load_adj_matrix();
    clock_t start = clock();
    RandomizedAlgorithm();
    //printf("\nk_hat: %d\n", k_hat);

    //printf("Inital Edges: %d mHat: %d Total edges: %d edges Removed: %d", initialEdges, m_hat, m_hat + cliqueEdges, edgesRemoved);
    int totalEdges = m_hat + cliqueEdges;
    float compression_ratio = compressionRatio();
    //printf("\nCompression Ratio: %f", r);
    clock_t stop = clock();
    double elapsed = ((double)(stop - start)) / CLOCKS_PER_SEC * 1000.0;
    //printf("\nExecution time: %lf", elapsed);
    printf("\n%d,%d,%d,%lf, %lf,%s\n", nodes, density, exp, compression_ratio, elapsed, cores);
    freeSet(S);
    free(S);
    getDeAllocate(graphNodes, adj_matrix);
    //getDeAllocateDouble(graphNodes, similarityMagnititude);
    free(leftCliqueSize);
    free(rightCliqueSize);
	free(S);
    //fclose(similarityArray);
    return 0;
}

