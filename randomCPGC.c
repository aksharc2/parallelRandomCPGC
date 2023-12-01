#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <stdbool.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <string.h>



void createFolder(const char* folderName, int rank) {
    // Check if the folder exists, create it if not
    struct stat st;
    if (stat(folderName, &st) == -1) {
        if (rank == 0) {
            if (mkdir(folderName, 0777) != 0) {
                perror("Error creating the folder");
                MPI_Abort(MPI_COMM_WORLD, 1); // Terminate MPI
            }
        }
    }
}



int get_k_hat(int graph_nodes, int m_hat, int delta) {
	int k_hat;
    float de = (2 * pow((double)graph_nodes, 2)) / m_hat;
    float nu = delta * log2((double)graph_nodes);
    k_hat = floor((double)nu / (log2((double)de)));
	return k_hat;
}

int main(int argc, char* argv[]) {
	
    MPI_Init(&argc, &argv);
	double start = MPI_Wtime();
    
    int comm_size;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    
	int delta = 1;
	int graphNodes = 8192 ; //atoi(argv[2]);
	const char* folderName = "temporaryEdges"; // Specify the folder name
    int arraySize = (int) ceil((double)(graphNodes * graphNodes)/(double)comm_size);
    // Get my rank
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	const char* filename = argv[1];
	
    // Create the window for a 10x10 matrix
    const int ROWS = graphNodes;
    const int COLS = graphNodes/comm_size;
    bool* adjMatrix_buffer;
    MPI_Win adjMatrix_window;
    MPI_Win_allocate_shared(arraySize * sizeof(bool), sizeof(bool), MPI_INFO_NULL, MPI_COMM_WORLD, &adjMatrix_buffer, &adjMatrix_window);
	MPI_Barrier(MPI_COMM_WORLD);
	
	// Initilizing the adjMatrix
	if (my_rank == 0){
		int lp, mp , rp, edges;
		int m_hat = 0;
		
		int adjOffset = 1;
		FILE *file = fopen(filename, "r");
		printf("started reading. \n");
		// Check if the file could be opened
		if (file == NULL) {
			printf("Failed to open the file.\n");
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		// Read the header information
		char line[256];
		fgets(line, sizeof(line), file);
		if (strncmp(line, "%%MatrixMarket matrix coordinate", 31) != 0) {
			printf("Invalid Matrix Market file.\n");
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		fgets(line, sizeof(line), file);
		while (line[0] == '%') {
			fgets(line, sizeof(line), file);
			//printf("%S", line);
		}

		// Parse the matrix size and number of non-zero values
		sscanf(line, "%d %d %d %d \n", &lp, &mp , &rp, &edges);
		printf("%d %d %d \n", lp, rp, edges);
			
		for (int i = 0; i < edges; i++) {
			int row, col;
			fscanf(file, "%d %d \n", &row, &col); // , &temp
			adjMatrix_buffer[ ((col - adjOffset) * graphNodes) + (row - adjOffset) ] = 1;
			m_hat++;
		}
		printf("Edges: %d\n",m_hat);
		fclose(file);
	// }
        
	// if (my_rank == 0) {
	    // printf("The array size is %d.\n", arraySize);
		
		
		// Print the initial matrix
		// for(int j = 0; j < comm_size * arraySize; j++){
			// printf("%d ", adjMatrix_buffer[j]);
			// if ((j+1)%graphNodes == 0)
				// printf("\n");
		// }
		// printf("\n");
		// if (lp >= rp)
			// graphNodes = lp;
		// else
			// graphNodes = rp;
		printf("Total edges : %d\n", m_hat);
		
// Algorithm Implementation starts here

		int S[graphNodes]; // initialized the set S
		for (int i  = 0; i < graphNodes; i++)
			S[i] = i;
		int setSize = graphNodes;
		createFolder(folderName, my_rank);
		char f_name[100];
		sprintf(f_name, "cliqueEdgesProcessor_%d.mtx", my_rank);
		const char* fileName = f_name; // = "output.txt";  // Specify the file name
		int *edgesRemoved = (int*)malloc((comm_size) * sizeof(int));
		MPI_Request gRequest, sRequest;		
		int k_hat = get_k_hat(graphNodes, m_hat, delta);
		printf("k_hat: %d\n", k_hat);
		for ( int r = 1; r < comm_size; r++){
			MPI_Isend(&k_hat, 1, MPI_INT, r, 0, MPI_COMM_WORLD, &sRequest);
		}
		int send_buf[k_hat + 2]; // = (int*)malloc((k_hat + 1) * sizeof(int));
		int q = graphNodes;
		while (k_hat > 1) {
			// MPI_Bcast(&k_hat, 1, MPI_INT, 0, MPI_COMM_WORLD);
			send_buf[0] = k_hat;
			int R = comm_size;
			int sendSize = k_hat + 2;
			while( setSize > 0){
				int edgeRemoved = 0;
				for ( int r = 1; r < comm_size; r++){
					send_buf[1] = q;
					q++;
					for (int k = 2; k <= k_hat + 1; k++){
						int randomVertex = rand() % setSize;
						send_buf[k]	= S[randomVertex];			
						for (int j = randomVertex; j < setSize - 1; j++){
							S[j] = S[j + 1];
						}
						setSize--;
						printf("				MPI process %d sends value %d to processor %d.\n", my_rank, send_buf[k], r);
					}	
					
					MPI_Send(&send_buf, sendSize, MPI_INT, r, 0, MPI_COMM_WORLD); // , &request			 , &sRequest		
					if (setSize <= k_hat){
						// printf("				MPI process %d finds common neighbour for graphNodes - s vertices.\n", my_rank);
						printf("**************setSize: %d\n", setSize);
						R = r + 1;
						int* targetIdx = (int*)malloc(setSize * sizeof(int));
						int neighboursIdx = 0;
						int* neighbours = (int*)malloc(graphNodes * sizeof(int));
						int remainingVertices = setSize;
						for(int k = 0; k < remainingVertices; k++){
							int p = (int) floor((S[k] * graphNodes)/arraySize);
							int idx = (S[k] * graphNodes) % arraySize;
							targetIdx[k] = (p - my_rank) * arraySize + idx;
							setSize--;
							printf("**************setSize: %d %d \n", setSize, k);
						}
						printf("**************setSize: %d\n", setSize);
						for(int i = 0; i < graphNodes; i++){
							int temp = 1;
							for(int k = 0; k < remainingVertices; k++){
								temp = temp * adjMatrix_buffer[targetIdx[k] + i];
							}
							if (temp){
								neighbours[neighboursIdx] = i;	
								// printf("[MPI process %d]----------- vertex u_%d has a common edge with (", my_rank, i);
									// for(int v = 0; v < remainingVertices; v++)
										// printf("v_%d, ", recv_buf[v]);
								// printf(")\n");
								edgeRemoved += remainingVertices;
								neighboursIdx++;
							}
						}					
						
						char filePath[256]; // Adjust the size as needed
						snprintf(filePath, sizeof(filePath), "%s/%s", folderName, fileName);
						FILE* file = NULL;
						file = fopen(filePath, "a");
						
						
						for(int i = 0; i < neighboursIdx; i++)
							fprintf(file, "%d %d\n", neighbours[i], q);
						
						
						for(int v = 0; v < remainingVertices; v++)
							fprintf(file, "%d %d\n", q, S[v]);			
						
						fclose(file);
			
						free(targetIdx);
						free(neighbours);
						q++;
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
				printf("remaining edges : %d\n", m_hat);
				
			}
			// calculate k_hat
			printf("remaining edges : %d\n", m_hat);
			k_hat = get_k_hat(graphNodes, m_hat, delta);	
			printf("k_hat: %d\n", k_hat);
		}
		
		printf("Sending terminating request\n");
		send_buf[0] = k_hat;

		for ( int r = 1; r < comm_size; r++){
			MPI_Isend(&send_buf, k_hat, MPI_INT, r, 0, MPI_COMM_WORLD, &sRequest); // , &request					
		}

		MPI_Wait(&sRequest, MPI_STATUS_IGNORE);

	}
	else{
		MPI_Request request;
		int received;
		char f_name[100];
		sprintf(f_name, "cliqueEdges_Processor_%d.mtx", my_rank);
		const char* fileName = f_name; // = "output.txt";  // Specify the file name
		int k_hat_buffer;
		MPI_Recv(&k_hat_buffer, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		int recv_buf[k_hat_buffer + 2]; // = (int*)malloc((k_hat + 1) * sizeof(int));
		while(1){
			int edgeRemoved = 0;
			// MPI_Bcast(&k_hat, 1, MPI_INT, 0, MPI_COMM_WORLD);				
			MPI_Recv(&recv_buf, k_hat_buffer + 2, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			int k_hat = recv_buf[0];
			int Q = recv_buf[1];
			printf("MPI process %d finding common neighbours for k_hat %d vertices.\n", my_rank, k_hat);
			if (k_hat < 2)
				break;
			
			printf("[MPI process %d]----------- received vertex (", my_rank);
			for(int v = 0; v <= k_hat + 1; v++)
				printf("%d, ", recv_buf[v]);
			printf(")\n");
			int* targetIdx = (int*)malloc(k_hat * sizeof(int));	
			int* neighbours = (int*)malloc(graphNodes * sizeof(int));				
			for(int k = 0; k < k_hat + 2; k++){
				received = recv_buf[k + 2];
				int p = (int) floor((received * graphNodes)/arraySize);
				int idx = (received * graphNodes) % arraySize;
				targetIdx[k] = (p - my_rank) * arraySize + idx;
			}
			int neighboursIdx = 0;
			for(int i = 0; i < graphNodes; i++){
				int temp = 1;
				for(int k = 0; k < k_hat; k++){
					temp = temp * adjMatrix_buffer[targetIdx[k] + i];
				}
				if (temp){
					neighbours[neighboursIdx] = i;					
					printf("[MPI process %d]----------- vertex u_%d has a common edge with (", my_rank, neighbours[neighboursIdx]);
						for(int v = 2; v <= k_hat+1; v++)
							printf("v_%d, ", recv_buf[v]);
					printf(")\n");
					neighboursIdx++;
					edgeRemoved += k_hat;
				}
			}
			MPI_Isend(&edgeRemoved, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &request);
			
			
			char filePath[256]; // Adjust the size as needed
			snprintf(filePath, sizeof(filePath), "%s/%s", folderName, fileName);
			FILE* file = NULL;
			file = fopen(filePath, "a");
			
			
			for(int i = 0; i < neighboursIdx; i++)
				fprintf(file, "%d %d\n", neighbours[i], Q);
			
			
			for(int v = 2; v <= k_hat+1; v++)
				fprintf(file, "%d %d\n", Q, recv_buf[v]);			
			
			fclose(file);
			
			
			
			free(targetIdx);
			free(neighbours);
		}
		// free(recv_buf);
		printf("Processor %d terminating as (k_hat) < 2 \n", my_rank);
	}

    MPI_Barrier(MPI_COMM_WORLD);
    // Destroy the window
    printf("[MPI process %d] adjMatrix_window destroyed.\n", my_rank);
    MPI_Win_free(&adjMatrix_window);
	// MPI_Win_free(&whileTerminator);
	if (my_rank == 0){
		double end = MPI_Wtime();
		printf("Time elapsed during the job: %.2fs.\n", end - start);
	}
    MPI_Finalize();

    return EXIT_SUCCESS;
}

