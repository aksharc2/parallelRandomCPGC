#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <ctype.h>
#include <string.h>
#include <inttypes.h> 
#include <time.h>

#define MAX_VERTICES 100000

bool isCompressed = false;
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
    graph->adjacencyList = (Node**)malloc(numVertices * sizeof(Node*));
    if (!graph->adjacencyList) {
        free(graph);
        return NULL;
    }

    for (uint32_t i = 0; i < numVertices; i++) {
        graph->adjacencyList[i] = NULL;
    }
    return graph;
}

void addEdge(Graph* graph, uint32_t src, uint32_t dest, bool isDirected) {
    Node* newNode = (Node*)malloc(sizeof(Node));
    if (!newNode) return;

    newNode->vertex = dest;
    newNode->next = graph->adjacencyList[src];
    graph->adjacencyList[src] = newNode;

    if (!isDirected) {
        Node* newNodeReverse = (Node*)malloc(sizeof(Node));
        if (!newNodeReverse) return;

        newNodeReverse->vertex = src;
        newNodeReverse->next = graph->adjacencyList[dest];
        graph->adjacencyList[dest] = newNodeReverse;
    }
}

// Queue implementation
typedef struct Queue {
    uint32_t* items;
    size_t front;
    size_t rear;
    size_t capacity;
} Queue;

Queue* createQueue(size_t capacity) {
    Queue* queue = (Queue*)malloc(sizeof(Queue));
    if (!queue) return NULL;

    queue->items = (uint32_t*)malloc(capacity * sizeof(uint32_t));
    if (!queue->items) {
        free(queue);
        return NULL;
    }
    queue->front = 0;
    queue->rear = 0;
    queue->capacity = capacity;
    return queue;
}

bool isEmpty(Queue* queue) {
    return queue->front == queue->rear;
}

void enqueue(Queue* queue, uint32_t value) {
    if ((queue->rear + 1) % queue->capacity == queue->front) {
        printf("Queue is full\n");
        return;
    }
    queue->items[queue->rear] = value;
    queue->rear = (queue->rear + 1) % queue->capacity;
}

uint32_t dequeue(Queue* queue) {
    if (isEmpty(queue)) {
        printf("Queue is empty\n");
        return 0xFFFF;
    }
    uint32_t item = queue->items[queue->front];
    queue->front = (queue->front + 1) % queue->capacity;
    return item;
}

void allPairsShortestPathsBFS(Graph* graph) {
    
    uint32_t numVertices =  graph->numVertices;
    int** distances = (int**)malloc(numVertices * sizeof(int*));
    if(!distances){
        fprintf(stderr, "Error: Memory allocation failed for distances array.\n");
        return;
    }
    for(uint32_t i = 0; i < numVertices; i++){
        distances[i] = (int*)malloc(numVertices * sizeof(int));
        if(!distances[i]){
            fprintf(stderr, "Error: Memory allocation failed for distances array.\n");
            for(uint32_t j = 0; j < i; j++){
                free(distances[j]);
            }
            free(distances);
            return;
        }
        for (int j = 0; j < numVertices; j++){
            distances[i][j] = -numVertices;
        }
    }

    for (uint32_t startVertex = 0; startVertex < numVertices; startVertex++) {
        bool* visited = (bool*)malloc(numVertices * sizeof(bool)); // Allocate on heap
        if (!visited) {
            fprintf(stderr, "Error: Memory allocation failed for visited array.\n");
            return; // Or handle the error appropriately
        }
        for (uint32_t i = 0; i < numVertices; i++){
            visited[i] = false;
        }

        Queue* queue = createQueue(numVertices + 1);
        if (!queue) {
            free(visited);
            return;
        }

        visited[startVertex] = true;
        enqueue(queue, startVertex);

        distances[startVertex][startVertex] = 0;
        int level = 0;

        while (!isEmpty(queue)) {
            size_t levelSize = queue->rear - queue->front;
            if (queue->rear < queue->front) {
                levelSize = queue->capacity - queue->front + queue->rear;
            }
            level++;
            for (size_t i = 0; i < levelSize; i++) {
                uint32_t currentVertex = dequeue(queue);

                Node* current = graph->adjacencyList[currentVertex];
                while (current != NULL) {
                    if (!visited[current->vertex]) {
                        visited[current->vertex] = true;
                        enqueue(queue, current->vertex);
                        distances[startVertex][current->vertex] = level;
                    }
                    current = current->next;
                }
            }
        }
        free(queue->items);
        free(queue);
        free(visited); // Free the allocated memory
    }

    
    for (uint32_t i = 0; i < numVertices; i++) {
        for (uint32_t j = 0; j < numVertices; j++) {
            totalDistance += distances[i][j];
            // printf("%d ", distances[i][j]);
            if (distances[i][j] > 0 && i < V && j < V){
                totalConnected += 1;
            }
        }
        // printf("\n");
    }
}

void freeGraph(Graph* graph) {
    for (uint32_t i = 0; i < graph->numVertices; i++) {
        Node* current = graph->adjacencyList[i];
        while (current != NULL) {
            Node* temp = current;
            current = current->next;
            free(temp);
        }
    }
    free(graph->adjacencyList);
    free(graph);
}

Graph* readMTXFile(const char* filename) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        perror("Unable to open file");
        exit(EXIT_FAILURE);
    }

    char line[256];
    while (fgets(line, sizeof(line), file)) {
        if (line[0] == '%' || line[0] == '\n') continue;
        break;
    }

    uint32_t edges, V_temp;
    uint32_t offset = 0;
    if (isCompressed){
        // sscanf(line, "%lu %lu %lu %lu", &V, &cliques, &V_temp, &edges);
        sscanf(line, "%" SCNu32 " %" SCNu32 " %" SCNu32 " %" SCNu32, &V, &cliques, &V_temp, &edges);
        offset = V + cliques;  // cliques is initialized to 0. So, for uncompressed graph it just adds 0.
    }
    else{
        // sscanf(line, "%lu %lu %lu", &V, &V_temp, &edges);
        sscanf(line, "%" SCNu32 " %" SCNu32 " %" SCNu32, &V, &V_temp, &edges);
    }

    if (V < V_temp){
        V = V_temp;
    }

    given_vertices = V;
    total_vetrices =  V + cliques;  // cliques is initialized to 0. So, for uncompressed graph it just adds 0.
    
    Graph* graph = createGraph(total_vetrices);

    while (fgets(line, sizeof(line), file)) {
        uint32_t src, dest, dest_prime;
        sscanf(line, "%u %u", &src, &dest_prime);
        src--; dest_prime--; // Convert to 0-based index

        if (dest_prime >= V){
            dest = dest_prime;
        }
        else{
            dest = dest_prime; // offset + 
        }
        addEdge(graph, src, dest, isDirected);
    }
    fclose(file);
    return graph;
}

int main(int argc, char* argv[]) {
    if (argc < 2) { // Minimum: program name and filename
        printf("Usage: %s <filename> [-d <0|1>] [-c <0|1>]\n", argv[0]);
        printf("  -d 0: Undirected, -d 1: Directed\n");
        printf("  -c 0: Uncompressed, -c 1: Compressed\n");
        return 1;
    }

    const char* filename = argv[1]; // Filename should always be the second argument
    bool d_found = false;
    bool c_found = false;

    for (int i = 2; i < argc; i++) {
        if (strcmp(argv[i], "-d") == 0) {
            if (i + 1 < argc && strlen(argv[i + 1]) == 1 && isdigit(argv[i + 1][0])) {
                isDirected = (argv[i + 1][0] == '1');
                i++; // Consume the value
                d_found = true;
            } else {
                fprintf(stderr, "Error: Invalid value for -d. Use 0 or 1.\n");
                return 1;
            }
        } else if (strcmp(argv[i], "-c") == 0) {
            if (i + 1 < argc && strlen(argv[i + 1]) == 1 && isdigit(argv[i + 1][0])) {
                isCompressed = (argv[i + 1][0] == '1');
                i++; // Consume the value
                c_found = true;
            } else {
                fprintf(stderr, "Error: Invalid value for -c. Use 0 or 1.\n");
                return 1;
            }
        } else {
            fprintf(stderr, "Error: Unknown argument: %s\n", argv[i]);
            return 1;
        }
    }

    // if (!d_found){
    //     fprintf(stderr, "Warning: -d argument not provided. Defaulting to Undirected\n");
    // }

    // if (!c_found){
    //     fprintf(stderr, "Warning: -c argument not provided. Defaulting to Uncompressed\n");
    // }

    clock_t start_time_readMTXFile = clock();
    Graph* graph = readMTXFile(argv[1]);
    double execution_time_readMTXFile = (double)(clock() - start_time_readMTXFile) / CLOCKS_PER_SEC;
    // printf("File read and graph structure execution time: %f seconds\n", execution_time_readMTXFile);

    clock_t start_time = clock();
    allPairsShortestPathsBFS(graph);
    double execution_time = (double)(clock() - start_time) / CLOCKS_PER_SEC;
    // printf("APSP execution time: %f seconds\n", execution_time);
    // printf("Total execution time: %f seconds\n", execution_time_readMTXFile + execution_time);
    // printf("All-Pairs Shortest Distances: %u\n", totalDistance);

    printf("%f, %f, %f, %u, %u\n", execution_time_readMTXFile, execution_time, execution_time_readMTXFile + execution_time, totalDistance, totalConnected);
    freeGraph(graph);

    return 0;
}