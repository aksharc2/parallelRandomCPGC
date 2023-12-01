#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int comm_size;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    int graphNodes = 32;
	int k_hat = 3;
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
            adjMatrix_buffer[i] = rand() % 2;// Initialize with unique values        
    }
    // Modify elements of the shared matrix
    MPI_Barrier(MPI_COMM_WORLD);
        
	if (my_rank == 0) {
	    printf("The array size is %d.\n", arraySize);
		// Print the initial matrix
		int m_hat = 0;
		for(int j = 0; j < comm_size * arraySize; j++){
			printf("%d ", adjMatrix_buffer[j]);
			m_hat = m_hat + adjMatrix_buffer[j];
			if ((j+1)%graphNodes == 0)
				printf("\n");
		}
		printf("\n");
		printf("Total edges : %d\n", m_hat);
		int *edgesRemoved = (int*)malloc((comm_size) * sizeof(int));
		
		MPI_Request gRequest, sRequest;
		
		
		int sendSize = k_hat +1 ;
		int send_buf[5];// = (int*)malloc((k_hat + 1) * sizeof(int));
		while (k_hat > 1) {
			// MPI_Bcast(&k_hat, 1, MPI_INT, 0, MPI_COMM_WORLD);
			int s = 0 ;
			int edgeRemoved = 0;
			send_buf[0] = k_hat;
			int R = comm_size;
			while( s < graphNodes){
				for ( int r = 1; r < comm_size; r++){
					for (int k = 1; k <= k_hat; k++){
						// int randomVertex = rand() % graphNodes;	
						send_buf[k] = s;
						printf("				MPI process %d sends value %d to processor %d.\n", my_rank, s, r);
						s++;
					}					
					MPI_Send(&send_buf, sendSize, MPI_INT, r, 0, MPI_COMM_WORLD); // , &request			 , &sRequest		
					if (s + k_hat > graphNodes && s != graphNodes){
						printf("				MPI process %d finds common neighbour for graphNodes - s vertices.\n", my_rank);
						R = r + 1;
						int vertices = graphNodes - s;
						int v = s;
						s = graphNodes;
						int* targetIdx = (int*)malloc(vertices * sizeof(int));
						int* recv_buf = (int*)malloc(vertices * sizeof(int));						
						for(int k = 0; k < vertices; k++){
							recv_buf[k] = v;
							int p = (int) floor((v * graphNodes)/arraySize);
							int idx = (v * graphNodes) % arraySize;
							targetIdx[k] = (p - my_rank) * arraySize + idx;
							v++;
						}
						for(int i = 0; i < graphNodes; i++){
							int temp = 1;
							for(int k = 0; k < vertices; k++){
								temp = temp * adjMatrix_buffer[targetIdx[k] + i];
							}
							if (temp){
								printf("[MPI process %d]----------- vertex u_%d has a common edge with (", my_rank, i);
									for(int v = 0; v < vertices; v++)
										printf("v_%d, ", recv_buf[v]);
								printf(")\n");
								edgeRemoved++;
							}
						}
						free(targetIdx);
						break;
					}
					// MPI_Wait(&sRequest, MPI_STATUS_IGNORE);
				}
				// MPI_Gather(&edgeRemoved, 1, MPI_INT, edgesRemoved, 1, MPI_INT, 0, MPI_COMM_WORLD); //, &gRequest
			    // MPI_Wait(&gRequest, MPI_STATUS_IGNORE);
				edgesRemoved[0] = edgeRemoved;
				m_hat = m_hat - edgesRemoved[0];
				printf("Edges removed by processor %d: %d \n", 0, edgesRemoved[0]);
				for ( int r = 1; r < R; r++){
					MPI_Recv(&edgeRemoved, 1, MPI_INT, r, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					edgesRemoved[r] = edgeRemoved;
					m_hat = m_hat - edgesRemoved[r];
					printf("Edges removed by processor %d: %d \n", r, edgesRemoved[r]);
				}
				printf("\n");
				printf("Total edges : %d\n", m_hat);
			}
			// calculate k_hat
			printf("Total edges : %d\n", m_hat);
			k_hat = 1;
			
		}
		printf("Sending terminating request\n");
		send_buf[0] = k_hat;
		for ( int r = 1; r < comm_size; r++){
			MPI_Isend(&send_buf, k_hat, MPI_INT, r, 0, MPI_COMM_WORLD, &sRequest); // , &request					
		}
		MPI_Wait(&sRequest, MPI_STATUS_IGNORE);
		// free(send_buf);
		// MPI_Bcast(&k_hat, 1, MPI_INT, 0, MPI_COMM_WORLD);
	}
	else{
		MPI_Request request;
		int received;
		int recv_buf[5]; // = (int*)malloc((k_hat + 1) * sizeof(int));
		while(1){
			int edgeRemoved = 0;
			// MPI_Bcast(&k_hat, 1, MPI_INT, 0, MPI_COMM_WORLD);				
			MPI_Recv(&recv_buf, (k_hat+1), MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			k_hat = recv_buf[0];
			printf("MPI process %d waiting with k_hat %d.\n", my_rank, k_hat);
			if (k_hat < 2)
				break;
			int* targetIdx = (int*)malloc(k_hat * sizeof(int));			
			for(int k = 0; k < k_hat; k++){
				received = recv_buf[k + 1];
				int p = (int) floor((received * graphNodes)/arraySize);
				int idx = (received * graphNodes) % arraySize;
				targetIdx[k] = (p - my_rank) * arraySize + idx;
			}
			for(int i = 0; i < graphNodes; i++){
				int temp = 1;
				for(int k = 0; k < k_hat; k++){
					temp = temp * adjMatrix_buffer[targetIdx[k] + i];
				}
				if (temp){
					printf("[MPI process %d]----------- vertex u_%d has a common edge with (", my_rank, i);
						for(int v = 1; v <= k_hat; v++)
							printf("v_%d, ", recv_buf[v]);
					printf(")\n");
					edgeRemoved++;
				}
			}
			MPI_Isend(&edgeRemoved, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &request);
			free(targetIdx);
		}
		// free(recv_buf);
		printf("Processor %d terminating as (k_hat) %d < 2 \n", my_rank, k_hat);
	}

    MPI_Barrier(MPI_COMM_WORLD);
    // Destroy the window
    printf("[MPI process %d] adjMatrix_window destroyed.\n", my_rank);
    MPI_Win_free(&adjMatrix_window);
	// MPI_Win_free(&whileTerminator);
    MPI_Finalize();

    return EXIT_SUCCESS;
}

