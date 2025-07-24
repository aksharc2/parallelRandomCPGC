/*On windows follow this instruction to get a running time:
  1.clock_t start = clock();
  2.clock_t stop = clock();
  3.elapsed = ((double)(stop - start)) / CLOCKS_PER_SEC * 1000.0;
 */

 /*On Grid follow this instruction to get a running time:
   0.struct timespec begin, end;
   1.clock_gettime(CLOCK_REALTIME, &begin);
   2.clock_gettime(CLOCK_REALTIME, &end);         
   3.long seconds = end.tv_sec - begin.tv_sec;
   4.long nanoseconds = end.tv_nsec - begin.tv_nsec;
   5.elapsed = (seconds + nanoseconds * 1e-9) * 1000;
  */

#pragma warning(disable:4996)
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <string.h>
#include <dirent.h>
#include <time.h>


#define MAXCHAR 256

void createFolder(const char* folderName) {
    // Check if the folder exists, create it if not
    struct stat st;
    if (stat(folderName, &st) == -1) {
		mkdir(folderName, 0777);
    }
}

struct timespec begin, end, ioTimeStart, ioTimeEnd, ioElapsed;
double elapsed_read;
double elapsed_write = 0.000;
bool** adj_matrix;
int* leftClique;
int* rightClique;
int k_hat;
float delta;
int m_hat = 0;
int initialEdges;
int edgesRemoved = 0;
int cliqueEdges = 0;
int q; // Clique index
int* S;
FILE* compressedFile;
int graphNodes;
struct timespec begin, end;

void load_adj_matrix(const char* inputFilename ){
	clock_gettime(CLOCK_REALTIME, &ioTimeStart);
    FILE* filePointer = fopen(inputFilename, "r");
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
			// printf("A(%d, %d) = %d \n", u, w, adj_matrix[u-1][w-1]);
		}
		initialEdges = m_hat = edges;
        fclose(filePointer);
    }
	clock_gettime(CLOCK_REALTIME, &ioTimeEnd);
	long seconds = ioTimeEnd.tv_sec - ioTimeStart.tv_sec;
    long nanoseconds = ioTimeEnd.tv_nsec - ioTimeStart.tv_nsec;
    elapsed_read = (seconds + nanoseconds * 1e-9) * 1000;
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

void saveCliqueEdges(int commonNeighbours){
	clock_gettime(CLOCK_REALTIME, &ioTimeStart);
	for(int u = 0; u < commonNeighbours; u++){
		fprintf(compressedFile, "%d %d\n", leftClique[u], q);
	}
	for (int w = 0; w < k_hat; w++){
		fprintf(compressedFile, "%d %d\n", q, rightClique[w]);
	}
	clock_gettime(CLOCK_REALTIME, &ioTimeEnd);
	long seconds = ioTimeEnd.tv_sec - ioTimeStart.tv_sec;
    long nanoseconds = ioTimeEnd.tv_nsec - ioTimeStart.tv_nsec;
    elapsed_write = elapsed_write + (seconds + nanoseconds * 1e-9) * 1000;
}

void findCommonNeighbours() {
    int i, j;
    // printf("\nLeft Partition: ");
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
            // printf("%d ", leftClique[u]);
            u++;
        }
    }
	if (u){
		cliqueEdges += u + k_hat;
		saveCliqueEdges(u);
	}
	//else{
		// printf("No common neighbours found\n");
		//for(int i = 0; i < k_hat; i++)
			//printf(" %d ", rightClique[i]);
	//}
}

float compressionRatio() {
    return ((float)initialEdges / (m_hat + cliqueEdges));
}

void RandomizedAlgorithm() {
    int nodeIndex;
	q = graphNodes;
	get_k_hat();
	// printf("k_hat: %d m_hat: %d cliqueEdges: %d", k_hat, m_hat, cliqueEdges);
    while (k_hat > 1) {
		for (int i = 0; i < graphNodes; i++){
			S[i] = i;
		}
		int setSize = graphNodes;
		while (setSize > 0){
			// printf("Clique %d \nright partition: ", q);
			leftClique = (int*)malloc(graphNodes * sizeof(int));
			rightClique = (int*)malloc(k_hat * sizeof(int));
			if (setSize >= k_hat){
				for(int i = 0; i < k_hat; i++){
					nodeIndex = rand() % setSize;
					rightClique[i] = S[nodeIndex];
					// printf("%d ", rightClique[i]);
					setSize--;
					for (int j = nodeIndex; j < setSize; j++){
						S[j] = S[j+1];
					}
				}
			}
			else{
				for(int i = 0; i < setSize; i++){
					rightClique[i] = S[i];
					// printf("%d ", rightClique[i]);
				}
				setSize = 0;
			}
			findCommonNeighbours();
			free(leftClique);
			free(rightClique);
			q++;
			// printf("k_hat: %d m_hat: %d cliqueEdges: %d\n", k_hat, m_hat, cliqueEdges);
		}
		get_k_hat();
    }
	// printf("k_hat: %d m_hat: %d cliqueEdges: %d\n", k_hat, m_hat, cliqueEdges);
}


int main(int argc, char* argv[]) {
	clock_gettime(CLOCK_REALTIME, &begin);
	const char* inputFilename = argv[1];
	graphNodes = atoi(argv[2]);
	delta = atof(argv[3]);
	int density = atoi(argv[4]);
	int instance = atoi(argv[5]);
    srand(time(NULL));
    adj_matrix = getAllocate(graphNodes); 
    S = (int*)malloc((graphNodes) * sizeof(int));
	char tempName[100];
	char compressedFileName[100];
	sprintf(tempName, "S_temp_%d_%d_%d.mtx", graphNodes, density, instance);
	sprintf(compressedFileName, "S_restrectured_graph_%d_%d_%d.mtx", graphNodes, density, instance);
	char filePath[256]; // Adjust the size as needed
	// snprintf(filePath, sizeof(filePath), "%s/%s", "sequentialCompressedGraphs", tempName);
    compressedFile = fopen(tempName, "w");
	
    load_adj_matrix(inputFilename);
    RandomizedAlgorithm();
	FILE * compFile = fopen(compressedFileName, "w");
	fclose(compressedFile);
    float edge_reduction_ratio = compressionRatio();
    //printf("\nCompression Ratio: %f", r);
	fprintf(compFile, "%%MatrixMarket matrix coordinate pattern general\n");
	fprintf(compFile, "%% Resulted restructured graph for given bipartite graph with %d nodes, %d density and %.1f delta.\n", graphNodes, density, delta);
	fprintf(compFile, "%d %d %d\n", graphNodes, graphNodes, m_hat + cliqueEdges);
	for(int i = 0; i < graphNodes; i++){
		for(int j = 0; j < graphNodes; j++){
			if (adj_matrix[i][j] != 0)
				fprintf(compFile, "%d %d\n", i, j);
		}
	}
	fclose(compFile);
	char command[250];
	sprintf(command, "cat %s >> %s", tempName, compressedFileName);
	system(command);
	getDeAllocate(graphNodes, adj_matrix);
	free(S);
    clock_gettime(CLOCK_REALTIME, &end);
	long seconds = end.tv_sec - begin.tv_sec;
    long nanoseconds = end.tv_nsec - begin.tv_nsec;
    double totalElapsed = (seconds + nanoseconds * 1e-9) * 1000;
	remove(tempName); 
    printf("%d,%d,%d, %.2f, %lf, %lf, %lf, %lf\n", graphNodes, density, instance, delta, edge_reduction_ratio, totalElapsed/1000, elapsed_read/1000, elapsed_write/1000);

    return 0;
}

