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
#include<stdio.h>
#include<string.h>
#include<stdbool.h>
#include<stdint.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>


#define MAXCHAR 100000


int** adj_matrix;
double** similarityMagnititude;
int** leftClique;
int** rightClique;
int* leftCliqueSize;
int* rightCliqueSize;
int k_hat;
int delta = 1;
int m_hat;
int initialEdges;
int edgesRemoved = 0;
int cliqueEdges = 0;
int q; // Clique index

int graphNodes;
struct timespec begin, end;
char f_name[100]; // name of the adjacency matrix file


FILE* similarityArray;


typedef struct {
    int* elements;
    int size;
    int capacity;
} Set;

Set* S;



void load_adj_matrix(){
    FILE* fpointer = fopen(f_name, "r");
    char line[MAXCHAR];
    if (!fpointer)
        printf("Cann't open the file\n");
    else {
        int row = 0;
        int col = 0;
        m_hat = 0;
        while (fgets(line, MAXCHAR, fpointer)) {
            col = 0;
            //strtok break down each line into smaller string
            char* value = strtok(line, ",");
            while (value != NULL) {
                adj_matrix[row][col] = atoi(value);
                m_hat += adj_matrix[row][col];
                value = strtok(NULL, ",");
                col++;
            }
            row++;
        }
        fclose(fpointer);
        initialEdges = m_hat;
    }
}



void get_k_hat() {
    float de = (2 * pow((double)graphNodes, 2)) / m_hat;
    float nu = delta * log2((double)graphNodes);
    k_hat = floor((double)nu / (log2((double)de)));
}


int** getAllocate(int n) {
    int i;
    int** arr = (int**)malloc((n) * sizeof(int*));
    for (i = 0; i < n; i++) {
        arr[i] = (int*)malloc((n) * sizeof(int));
    }
    return arr;
}

double** getAllocateDouble(int n) {
    int i;
    double** arr = (double**)malloc((n) * sizeof(double*));
    for (i = 0; i < n; i++) {
        arr[i] = (double*)malloc((n) * sizeof(double));
    }
    return arr;
}

void getDeAllocate(int n, int** arr) {
    int i;
    for (i = 0; i < n; i++) {
        free(arr[i]);
    }
    free(arr);
}

void getDeAllocateDouble(int n, double** arr) {
    int i;
    for (i = 0; i < n; i++) {
        free(arr[i]);
    }
    free(arr);
}





void initialize(Set* set, int nodes) {
    set->size = 0;
    set->capacity = nodes;  // Initial capacity
    set->elements = (int*)malloc(set->capacity * sizeof(int));
    if (set->elements == NULL) {
        printf("Memory allocation failed.\n");
        exit(EXIT_FAILURE);
    }
}

void insert(Set* set, int element) {
    // Check if the element already exists in the set
    //for (int i = 0; i < set->size; i++) {
    //    printf("Element %d in the set.\n", set->elements[i]);
    //    if (set->elements[i] == element) {
    //        printf("Element %d already exists in the set.\n", element);
    //        return;
    //    }
    //}

    // Resize the array if it's full
    //if (set->size == set->capacity) {
    //    set->capacity *= 2;  // Double the capacity
    //    set->elements = (int*)realloc(set->elements, set->capacity * sizeof(int));
    //    if (set->elements == NULL) {
    //        printf("Memory reallocation failed.\n");
    //        exit(EXIT_FAILURE);
    //    }
    //}

    // Insert the element into the set
    set->elements[set->size] = element;
    set->size++;
    //printf("Element %d inserted into the set.\n", element);
}

int removeElement(Set* set, int index) {
    //int found = 0;

    // Find the index of the element
    //int index;
    //for (index = 0; index < set->size; index++) {
    //    if (set->elements[index] == element) {
    //        found = 1;
    //        break;
    //    }
    //}

    //if (found) {
        // Shift elements to the left to fill the gap
    int element = set->elements[index];
    for (int i = index; i < set->size - 1; i++) {
        set->elements[i] = set->elements[i + 1];
    }

    set->size--;
    //printf("Element %d removed from the set.\n", element);
    return element;
    //}
    //else {
    //    printf("Element %d not found in the set.\n", element);
    //    return false;
    //}
}

void search(Set* set, int element) {
    for (int i = 0; i < set->size; i++) {
        if (set->elements[i] == element) {
            printf("Element %d found in the set.\n", element);
            return;
        }
    }
    printf("Element %d not found in the set.\n", element);
}

void display(Set* set) {
    printf("Elements in the set: ");
    for (int i = 0; i < set->size; i++) {
        printf("%d ", set->elements[i]);
    }
    printf("\n");
}

void freeSet(Set* set) {
    free(set->elements);
    set->size = 0;
    set->capacity = 0;
    set->elements = NULL;
}


void updateSet(Set* S) {
    int i;
    for (i = 0; i < graphNodes; i++)
        insert(S, i);
}


void findCommonNeighbours(int k) {
    int i, j;
    printf("\nLeft Partition: ");
    int u = 0;
    for (i = 0; i < graphNodes; i++) {
        int flag = 0;
        for (j = 0; j < rightCliqueSize[k] && !flag; j++) {
            if (adj_matrix[i][rightClique[k][j]] == 0) {
                flag = 1;
            }
        }
        if (!flag) {
            leftClique[k][u] = i;
            for (j = 0; j < rightCliqueSize[k] && !flag; j++) {
                adj_matrix[i][rightClique[k][j]] = 0;
                m_hat -= 1; // need to remove this for parallel implementation
                edgesRemoved += 1; // need to remove this for parallel implementation
            }
            printf("%d ", leftClique[k][u]);
            u++;
        }
    }
    leftCliqueSize[k] = u;
    if (u > 0)
        cliqueEdges += rightCliqueSize[k] + leftCliqueSize[k];
    //printf("\n Size of Left Partition of %d Clique: %d", d, u);
    //printf("\nEdges in Clique: %d", cliqueEdges);
}


float compressionRatio() {
    return ((float)initialEdges / (m_hat + cliqueEdges));
}




void RandomizedAlgorithm() {
    int nodeIndex;
    while (k_hat > 1) {
        int n = 0;
        printf("Clique %d \nright partition: ", q);
        while (S->size + n >= k_hat) { // this condition is for not picking the same vertex in the right clique of the last iteration
            nodeIndex = rand() % (S->size);
            rightClique[q][n] = removeElement(S, nodeIndex);
            printf(" %d ", rightClique[q][n]);
            if (n < k_hat - 1 )
                n++;
            else {
                
                rightCliqueSize[q] = n+1;
                findCommonNeighbours(q);
                q++;
                n = 0;
                printf("\n\nClique %d \nright partition: ", q);
            }
        }
        if (graphNodes % k_hat) {
            n = 0;
            while (S->size > 0) {
                rightClique[q][n] = removeElement(S, n);
                printf(" %d ", rightClique[q][n]);
                n++;
            }
            rightCliqueSize[q] = graphNodes % k_hat;
            findCommonNeighbours(q);
            q++;
        }
        get_k_hat();
        //printf("\nk_hat: %d\n", k_hat);
        if (S->size == 0 && k_hat > 1) {
            updateSet(S);
        }
        //printf("\nNew Set Formed\n");
        //display(S);
    }
}

// Cosine Similarity functions
double cosineSimilarity(int i, int j) {
    int size = graphNodes;
    int dot_product = 0;
    int magnitude1 = 0;
    int magnitude2 = 0;

    for (int n = 0; n < size; n++) {
        dot_product += adj_matrix[n][i] * adj_matrix[n][j];
        magnitude1 += adj_matrix[n][i];
        magnitude2 += adj_matrix[n][j];
    }

    double similarity = dot_product / (sqrt(magnitude1) * sqrt(magnitude2));
    return similarity;
}



//void getSimilarity() {
//    int i, j;
//    for (i = 0; i < graphNodes; i++) {
//        for (j = 0; j < graphNodes; j++) {
//            similarityMagnititude[i][j] = cosineSimilarity(i, j);
//            printf("%lf ", similarityMagnititude[i][j]);
//            fprintf(similarityArray, "%lf,", similarityMagnititude[i][j]);
//        }
//        printf("\n");
//        fprintf(similarityArray, "\n");
//    }
//}



int main() {
    char cores[8] = "randSeq";
    int nodes = 32; // atoi(argv[1]);  // int argc, char* argv[]
    int density = 80; // atoi(argv[2]);
    int exp = 1; // atoi(argv[3]);
    int multiplier;
    srand(time(NULL));
    graphNodes = nodes;
    multiplier = ceil(log10((double)graphNodes));
    S = (Set*)malloc(sizeof(Set));
    //graphNodes = 32;
    adj_matrix = getAllocate(graphNodes); 
    //similarityMagnititude = getAllocateDouble(graphNodes);
    
    rightClique = getAllocate(multiplier * graphNodes); //  Could lower the size of right partition array as at max we would extract K_hat vertices.
    leftClique = getAllocate(multiplier * graphNodes);
    leftCliqueSize = (int*)malloc((multiplier * graphNodes) * sizeof(int));
    rightCliqueSize = (int*)malloc((multiplier * graphNodes) * sizeof(int));

    //sprintf(f_name, "Bipartite%dX%d.csv", graphNodes, graphNodes);
    sprintf(f_name, "New_generated_data/Bipartite_%dX%d/%d/Bipartite_%dX%d_%d_%d.csv", nodes, nodes, density, nodes, nodes, density, exp);

    load_adj_matrix();
    clock_t start = clock();

    initialize(S, graphNodes);
    updateSet(S);

    get_k_hat();

    //similarityArray = fopen("SimilarityArray.csv", "w");
    
    q = 0;
    //display(S);
    //getSimilarity();
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
    getDeAllocate(graphNodes * multiplier, rightClique);
    getDeAllocate(graphNodes * multiplier, leftClique);
    //getDeAllocateDouble(graphNodes, similarityMagnititude);
    free(leftCliqueSize);
    free(rightCliqueSize);
    //fclose(similarityArray);
    return 0;
}

