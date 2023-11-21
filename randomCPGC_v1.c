#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    // Check that only 2 MPI processes are spawned
    int comm_size;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    int graphNodes = 5;
    int arraySize = (int) ceil((double)(graphNodes * graphNodes)/(double)comm_size);

    // Get my rank
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    // Create the window for a 10x10 matrix
    const int ROWS = graphNodes;
    const int COLS = graphNodes/comm_size;
    int* adjMatrix_buffer;
    MPI_Win adjMatrix_window, leftClique_window, rightClique_window, whileTerminator;
    MPI_Win_allocate_shared(arraySize * sizeof(int), sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &adjMatrix_buffer, &adjMatrix_window);
    printf("[MPI process %d] adjMatrix_window created.\n", my_rank);
    
   
    // Initialize the shared matrix
    for (int i = 0; i < arraySize; i++) {       
            adjMatrix_buffer[i] = 0 ; // Initialize with unique values        
    }


    // Print the initial matrix
    printf("[MPI process %d] Initial matrix:\n", my_rank);
	for (int i = 0; i < arraySize; i++) {  
	    printf("%d ", adjMatrix_buffer[i]);
	}
	printf("\n");
    // Modify elements of the shared matrix
    MPI_Barrier(MPI_COMM_WORLD);
        
	if (my_rank == 0) {
	    printf("The array size is %d.\n", arraySize);
	    MPI_Request request;
	    
	    for ( int r = 1; r < comm_size; r++){	      
		    int randomVertex = rand() % graphNodes;		    
		    printf("MPI process %d sends value %d to processor %d.\n", my_rank, randomVertex, r);
		    MPI_Isend(&randomVertex, 1, MPI_INT, r, 0, MPI_COMM_WORLD, &request);
            }
            // Do other things while the MPI_Isend completes
            // <...>
 
            // Let's wait for the MPI_Isend to complete before progressing further.
            //MPI_Wait(&request, MPI_STATUS_IGNORE);


	}
	else{
            int received;
            MPI_Recv(&received, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            printf("MPI process %d received value: %d.\n", my_rank, received);
            int p = (int) floor((received * graphNodes)/arraySize);
            int idx = (received * graphNodes) % arraySize;
			printf("Accesing processor %d and index %d \n", p, idx);
			int targetIdx = (p - my_rank) * arraySize + idx;
			printf("Processor %d target Id %d \n", my_rank, targetIdx);
			for (int i = targetIdx; i < (graphNodes +  targetIdx); i++) { 
				adjMatrix_buffer[i] = my_rank;
				//printf("%d ", adjMatrix_buffer[i]);
			}
	}

    MPI_Barrier(MPI_COMM_WORLD);

    // Print the updated matrix
    printf("[MPI process %d] Updated matrix:\n", my_rank);
 
    for (int i = 0; i < arraySize; i++) {  
    	printf("(%d) %d ", my_rank, adjMatrix_buffer[i]);
    }

    // Destroy the window
    printf("[MPI process %d] adjMatrix_window destroyed.\n", my_rank);
    MPI_Win_free(&adjMatrix_window);

    MPI_Finalize();

    return EXIT_SUCCESS;
}

