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
	int k_hat;
    int arraySize = (int) ceil((double)(graphNodes * graphNodes)/(double)comm_size);

    // Get my rank
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	
    // Create the window for a 10x10 matrix
    const int ROWS = graphNodes;
    const int COLS = graphNodes/comm_size;
	int termination = -99;
    int* adjMatrix_buffer;
	// int* whileTerminator_buffer;
    MPI_Win adjMatrix_window;
    MPI_Win_allocate_shared(arraySize * sizeof(int), sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &adjMatrix_buffer, &adjMatrix_window);
	// MPI_Win_allocate_shared(1 * sizeof(int), sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &whileTerminator_buffer, &whileTerminator);
   
	// *whileTerminator_buffer = 1;
    // Initialize the shared matrix
    for (int i = 0; i < arraySize; i++) {       
            adjMatrix_buffer[i] = 0 ; // Initialize with unique values        
    }
    // Modify elements of the shared matrix
    MPI_Barrier(MPI_COMM_WORLD);
        
	if (my_rank == 0) {
	    printf("The array size is %d.\n", arraySize);
		// Print the initial matrix
		
		for(int j = 0; j < comm_size; j++){
			for (int i = 0; i < arraySize; i++) {  
				printf("%d ", adjMatrix_buffer[j* comm_size + i]);
			}
			printf("\n");
		}
		
		MPI_Request request;
		int s = 0 ;
		k_hat = 2;
		MPI_Bcast(&k_hat, 1, MPI_INT, 0, MPI_COMM_WORLD);
	    while( s < graphNodes){
			for ( int r = 1; r < comm_size; r++){	      
				int randomVertex = rand() % graphNodes;		    
				printf("MPI process %d sends value %d to processor %d.\n", my_rank, randomVertex, r);
				MPI_Isend(&randomVertex, 1, MPI_INT, r, 0, MPI_COMM_WORLD, &request);
				s++;
				if (s == graphNodes)
					break;
			}
			//MPI_Barrier(MPI_COMM_WORLD);
		}
		k_hat = 1;
		MPI_Bcast(&k_hat, 1, MPI_INT, 0, MPI_COMM_WORLD);
		printf("Sending terminating request\n");
		// for (int c = 0; c < comm_size; c++){
			// MPI_Isend(&termination, 1, MPI_INT, c, 0, MPI_COMM_WORLD, &request);
		// }
		
	}
	else{
			MPI_Request request;
			int received;
			while(1){
				
				MPI_Bcast(&k_hat, 1, MPI_INT, 0, MPI_COMM_WORLD);
				printf("MPI process %d received k_hat value: %d.\n", my_rank, k_hat);
				if (k_hat < 2)
					break;
				MPI_Recv(&received, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				//MPI_Wait(&request, MPI_STATUS_IGNORE);

				printf("MPI process %d received value: %d.\n", my_rank, received);
				int p = (int) floor((received * graphNodes)/arraySize);
				int idx = (received * graphNodes) % arraySize;
				printf("Accesing processor %d and index %d \n", p, idx);
				int targetIdx = (p - my_rank) * arraySize + idx;
				printf("Processor %d target Id %d \n", my_rank, targetIdx);
				for (int i = targetIdx; i < (graphNodes +  targetIdx); i++) { 
					adjMatrix_buffer[i] = my_rank;
					printf("%d ", adjMatrix_buffer[i]);
				}
				printf("\n");
			}
			printf("Processor %d terminating as (k_hat) %d < 2 \n", my_rank, k_hat);
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
	// MPI_Win_free(&whileTerminator);
    MPI_Finalize();

    return EXIT_SUCCESS;
}

