#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <ctype.h>
#include <string.h>
#include <inttypes.h>
#include <time.h>

#define MAX_VERTICES 10000
#define INF 0x3f3f3f3f

bool isCompressed = true;
bool isDirected = false;
uint32_t cliques = 0;
uint32_t V;
uint32_t given_vertices;
uint32_t total_vetrices;
uint32_t totalDistance = 0;
uint32_t totalConnected = 0;

typedef struct Node {
    uint32_t vertex;
    struct Node* next;
} Node;

typedef struct Graph {
    Node** adjacencyList;
    uint32_t numVertices;
} Graph;

Graph* createGraph(uint32_t numVertices) {
    Graph* graph = (Graph*)malloc(sizeof(Graph));
    if (!graph) return NULL;

    graph->numVertices = numVertices;
    graph->adjacencyList = (Node**)calloc(numVertices, sizeof(Node*));
    if (!graph->adjacencyList) {
        free(graph);
        return NULL;
    }
    return graph;
}

void addEdge(Graph* graph, uint32_t src, uint32_t dest, bool isDirected) {
    if (src >= graph->numVertices || dest >= graph->numVertices) {
        fprintf(stderr, "Error: Invalid edge %u -> %u\n", src, dest);
        exit(EXIT_FAILURE);
    }

    Node* newNode = (Node*)malloc(sizeof(Node));
    newNode->vertex = dest;
    newNode->next = graph->adjacencyList[src];
    graph->adjacencyList[src] = newNode;

    if (!isDirected) {
        Node* newNodeReverse = (Node*)malloc(sizeof(Node));
        newNodeReverse->vertex = src;
        newNodeReverse->next = graph->adjacencyList[dest];
        graph->adjacencyList[dest] = newNodeReverse;
    }
}

// Queue implementation
typedef struct Queue {
    uint32_t* items;
    size_t front, rear, size, capacity;
} Queue;

Queue* createQueue(size_t capacity) {
    Queue* q = (Queue*)malloc(sizeof(Queue));
    q->items = (uint32_t*)malloc(capacity * sizeof(uint32_t));
    q->front = 0;
    q->rear = 0;
    q->size = 0;
    q->capacity = capacity;
    return q;
}

bool isEmpty(Queue* q) {
    return q->size == 0;
}

void enqueue(Queue* q, uint32_t value) {
    if (q->size == q->capacity) {
        fprintf(stderr, "Queue overflow\n");
        exit(EXIT_FAILURE);
    }
    q->items[q->rear] = value;
    q->rear = (q->rear + 1) % q->capacity;
    q->size++;
}

uint32_t dequeue(Queue* q) {
    if (isEmpty(q)) {
        fprintf(stderr, "Queue underflow\n");
        exit(EXIT_FAILURE);
    }
    uint32_t item = q->items[q->front];
    q->front = (q->front + 1) % q->capacity;
    q->size--;
    return item;
}

void allPairsShortestPathsBFS(Graph* graph) {
    uint32_t numVertices = graph->numVertices;

    int** distances = (int**)malloc(numVertices * sizeof(int*));
    for (uint32_t i = 0; i < numVertices; i++) {
        distances[i] = (int*)malloc(numVertices * sizeof(int));
        for (uint32_t j = 0; j < numVertices; j++)
            distances[i][j] = -1;
    }

    for (uint32_t start = 0; start < numVertices; start++) {
        bool* visited = (bool*)calloc(numVertices, sizeof(bool));
        Queue* q = createQueue(numVertices);

        visited[start] = true;
        distances[start][start] = 0;
        enqueue(q, start);

        while (!isEmpty(q)) {
            uint32_t u = dequeue(q);
            Node* curr = graph->adjacencyList[u];
            while (curr) {
                uint32_t v = curr->vertex;
                if (v >= numVertices) {
                    fprintf(stderr, "Error: Vertex %u out of bounds.\n", v);
                    exit(EXIT_FAILURE);
                }
                if (!visited[v]) {
                    visited[v] = true;
                    distances[start][v] = distances[start][u] + 1;
                    enqueue(q, v);
                }
                curr = curr->next;
            }
        }

        free(q->items);
        free(q);
        free(visited);
    }

    for (uint32_t i = 0; i < V; i++) {
        for (uint32_t j = 0; j < V; j++) {
            if (distances[i][j] > 0) {
                totalDistance += distances[i][j];
                totalConnected++;
            }
        }
    }

    for (uint32_t i = 0; i < numVertices; i++) {
        free(distances[i]);
    }
    free(distances);
}

void freeGraph(Graph* graph) {
    for (uint32_t i = 0; i < graph->numVertices; i++) {
        Node* curr = graph->adjacencyList[i];
        while (curr) {
            Node* tmp = curr;
            curr = curr->next;
            free(tmp);
        }
    }
    free(graph->adjacencyList);
    free(graph);
}

Graph* readMTXFile(const char* filename) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        perror("File open error");
        exit(EXIT_FAILURE);
    }

    char line[256];
    while (fgets(line, sizeof(line), file)) {
        if (line[0] != '%' && line[0] != '\n')
            break;
    }

    uint32_t edges, V_temp;
    if (isCompressed) {
        sscanf(line, "%" SCNu32 " %" SCNu32 " %" SCNu32 " %" SCNu32, &V, &cliques, &V_temp, &edges);
    } else {
        sscanf(line, "%" SCNu32 " %" SCNu32 " %" SCNu32, &V, &V_temp, &edges);
    }
    printf("%d, %d, %d, %d", V, cliques, V_temp, edges);
    if (V < V_temp) V = V_temp;
    given_vertices = V;
    total_vetrices = V + cliques + 1;

    Graph* graph = createGraph(total_vetrices);

    while (fgets(line, sizeof(line), file)) {
        uint32_t src, dest;
        if (sscanf(line, "%u %u", &src, &dest) != 2)
            continue;

        if (src >= total_vetrices || dest >= total_vetrices) {
            fprintf(stderr, "Invalid edge %u -> %u\n", src, dest);
            continue;
        }

        addEdge(graph, src, dest, isDirected);
    }

    fclose(file);
    return graph;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <filename> [-d <0|1>] [-c <0|1>]\n", argv[0]);
        return EXIT_FAILURE;
    }

    const char* filename = argv[1];

    for (int i = 2; i < argc; i++) {
        if (strcmp(argv[i], "-d") == 0 && i + 1 < argc) {
            isDirected = (argv[++i][0] == '1');
        } else if (strcmp(argv[i], "-c") == 0 && i + 1 < argc) {
            isCompressed = (argv[++i][0] == '1');
        } else {
            fprintf(stderr, "Unknown argument: %s\n", argv[i]);
            return EXIT_FAILURE;
        }
    }

    clock_t t1 = clock();
    Graph* graph = readMTXFile(filename);
    clock_t t2 = clock();
    double t_read = (double)(t2 - t1) / CLOCKS_PER_SEC;

    clock_t t3 = clock();
    allPairsShortestPathsBFS(graph);
    clock_t t4 = clock();
    double t_apsp = (double)(t4 - t3) / CLOCKS_PER_SEC;

    printf("\n%f, %f, %f, %u, %u\n", t_read, t_apsp, t_read + t_apsp, totalDistance, totalConnected);

    freeGraph(graph);
    return 0;
}
