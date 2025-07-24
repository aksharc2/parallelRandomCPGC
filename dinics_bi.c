/*This code is the implementation of Dinic's algorithm for bipartite graph*/
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include<time.h>

 

int* leftPartitionSize;
int middlePartitionSize = 0;
int* rightPartitionSize;
int* edges;
double run_time, total_run_time;
int maximumFlow;
 

char* filename[100];
int total_nodes;
struct timespec readingmtx, begin, end;

 

 

struct Edge {
    int v;
    int flow;
    int C;
    int rev;
};

 

struct Graph {
    int V;
    int* level;
    struct Edge** adj;
    int* edgeCount;
};

 

 

struct Graph* createGraph(int V) {
    struct Graph* graph = (struct Graph*)malloc(sizeof(struct Graph));
    graph->V = V;
    graph->adj = (struct Edge**)malloc(V * sizeof(struct Edge*));
    for (int i = 0; i < V; i++) {
        graph->adj[i] = (struct Edge*)malloc(V * sizeof(struct Edge));
    }
    graph->level = (int*)malloc(V * sizeof(int));
    graph->edgeCount = (int*)calloc(V, sizeof(int));
    return graph;
}

 

struct Graph* graph;

 

void addEdge(struct Graph* graph, int u, int v, int C) {
    graph->adj[u][graph->edgeCount[u]] = (struct Edge){ v, 0, C, graph->edgeCount[v] };
    graph->adj[v][graph->edgeCount[v]] = (struct Edge){ u, 0, 0, graph->edgeCount[u] };
    graph->edgeCount[u]++;
    graph->edgeCount[v]++;
}

 

bool BFS(struct Graph* graph, int s, int t) {
    for (int i = 0; i < graph->V; i++)
        graph->level[i] = -1;

 

    graph->level[s] = 0;
    int* q = (int*)malloc(graph->V * sizeof(int));
    int front = 0, rear = 0;
    q[rear++] = s;

 

    while (front < rear) {
        int u = q[front++];
        for (int i = 0; i < graph->edgeCount[u]; i++) {
            struct Edge e = graph->adj[u][i];
            if (graph->level[e.v] < 0 && e.flow < e.C) {
                graph->level[e.v] = graph->level[u] + 1;
                q[rear++] = e.v;
            }
        }
    }

 

    free(q);
    return graph->level[t] < 0 ? false : true;
}

 

int sendFlow(struct Graph* graph, int u, int flow, int t, int start[]) {
    if (u == t) {
        return flow;
    }

 

    for (; start[u] < graph->edgeCount[u]; start[u]++) {
        struct Edge* e = &graph->adj[u][start[u]];

 

        if (graph->level[e->v] == graph->level[u] + 1 && e->flow < e->C) {
            int curr_flow = flow < e->C - e->flow ? flow : e->C - e->flow;

 

            int temp_flow = sendFlow(graph, e->v, curr_flow, t, start);

 

            if (temp_flow > 0) {
                e->flow += temp_flow;
                graph->adj[e->v][e->rev].flow -= temp_flow;
                return temp_flow;
            }
        }
    }

 

    return 0;
}

 

int DinicMaxflow(struct Graph* graph, int s, int t) {
    if (s == t) {
        return -1;
    }

 

    int total = 0;

 

    while (BFS(graph, s, t) == true) {
        int* start = (int*)malloc(graph->V * sizeof(int));
        for (int i = 0; i < graph->V; i++)
            start[i] = 0;

 

        int flow;
        while ((flow = sendFlow(graph, s, INT_MAX, t, start)) > 0) {
            total += flow;
        }

 

        free(start);
    }

 

    return total;
}

 

 

void readMatrixMarketFile() {
    FILE* file = fopen(filename, "r");

 

    // Check if the file could be opened
    if (file == NULL) {
        printf("Failed to open the file.\n");
        exit(1);
    }

 

    // Read the header information
    char line[256];
    fgets(line, sizeof(line), file);
    if (strncmp(line, "%%MatrixMarket matrix coordinate", 31) != 0) {
        printf("Invalid Matrix Market file.\n");
        exit(1);
    }

 

    fgets(line, sizeof(line), file);
    while (line[0] == '%') {
        fgets(line, sizeof(line), file);
    }

 

    // Parse the matrix size and number of non-zero values
    sscanf(line, "%d %d %d\n", &leftPartitionSize, &rightPartitionSize, &edges);
    total_nodes = (int)leftPartitionSize + (int)rightPartitionSize + middlePartitionSize + 2;
    //printf("total_nodes: %d\n", total_nodes);
    graph = createGraph(total_nodes);
    int i;
    for (i = 0; i < edges; i++) {
        int row, col;
        fscanf(file, "%d %d \n", &row, &col);
        //addEdge(graph, row, col + (int)leftPartitionSize + (int)middlePartitionSize, 1);

 

 

        if (col > (int) rightPartitionSize)
            addEdge(graph, row, col, 1);
        else
            addEdge(graph, row, col + (int)leftPartitionSize + middlePartitionSize, 1);

 

 

        if (i < (int)leftPartitionSize) {
            //printf("%d %d \n", 0, i+1 );
            addEdge(graph, 0, i + 1, 1);
        }
        //printf("%d %d \n", row, col + (int)leftPartitionSize);

 

        if (i < (int)rightPartitionSize) {
            //printf("%d %d \n", i + 1 + (int)leftPartitionSize, total_nodes -1);
            addEdge(graph, i + 1 + (int)leftPartitionSize + middlePartitionSize, total_nodes - 1, 1);
        }
        //addEdge(graph, row, col, 1);
        //    //adjMatrix[row][col] = 1;
        //    //matrix.row_indices[i] = row;  // Convert to 0-based indexing
        //    //matrix.col_indices[i] = col;  // Convert to 0-based indexing
        //    //matrix.values[i] = 1;
    }

 

    fclose(file);
    //return matrix;
}



int main(int argc, char* argv[]) {
    int nodes;
    int density;
    int exp;
    float delta = 1.0;
    if (argc > 1) {
        nodes = atoi(argv[1]);
        density = atoi(argv[2]);
        exp = atoi(argv[3]);
	//delta = atof(argv[4]);

        sprintf(filename, "/ocean/projects/cis230093p/achavan/graphCompression/datasets/bipartite_graph_32768_%d_11.mtx", density);
        //sprintf(filename, "data3_1/bipartite_graph_%d_%d_%d.mtx", nodes, density, exp); bipartite_graph_32768_80_11.mtx
    }
    clock_gettime(CLOCK_REALTIME, &readingmtx);
    readMatrixMarketFile();
    clock_gettime(CLOCK_REALTIME, &begin);
    maximumFlow = DinicMaxflow(graph, 0, total_nodes - 1);
    //printf("Maximum flow: %d\n", DinicMaxflow(graph, 0, total_nodes - 1));
    clock_gettime(CLOCK_REALTIME, &end);
    long seconds = end.tv_sec - begin.tv_sec;
    long nanoseconds = end.tv_nsec - begin.tv_nsec;
    run_time = (seconds + nanoseconds * 1e-9) * 1000;
 

    long seconds2 = end.tv_sec - readingmtx.tv_sec;
    long nanoseconds2 = end.tv_nsec - readingmtx.tv_nsec;
    total_run_time = (seconds2 + nanoseconds2 * 1e-9) * 1000;

 
	
    printf("%d, %d, %d, %d, %d, %f, %f\n", nodes, total_nodes, density, exp, maximumFlow, run_time, total_run_time);

    //printf("Execution Time: %f ms\nTotal execution time: %f ms\n", run_time, total_run_time);

    // Free dynamically allocated memory
    for (int i = 0; i < graph->V; i++) {
        free(graph->adj[i]);
    }
    free(graph->adj);
    free(graph->level);
    free(graph->edgeCount);
    free(graph);

 

    return 0;
}
