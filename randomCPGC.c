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

int main(int argc, char* argv[]) {
	const char* folderName = "temporaryEdges"; // Specify the folder name
	createFolder(folderName);
	
    MPI_Init(&argc, &argv);
	double start = MPI_Wtime();
    
    int comm_size;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    
	float delta = atof(argv[3]);
	int graphNodes = atoi(argv[2]);

    int arraySize = (int) ceil((double)(graphNodes * graphNodes)/(double)comm_size);
    // Get my rank
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	const char* filename = argv[1];
	int edgesAddedToClique = 0;
	int TotalCliqueEdges = 0;
	int InitialEdges = 0;
	int m_hat = 0;
    // Create the window for a 10x10 matrix
    const int ROWS = graphNodes;
    const int COLS = graphNodes/comm_size;
    bool* adjMatrix_buffer;
    MPI_Win adjMatrix_window;
    MPI_Win_allocate_shared(arraySize * sizeof(bool), sizeof(bool), MPI_INFO_NULL, MPI_COMM_WORLD, &adjMatrix_buffer, &adjMatrix_window);
	MPI_Barrier(MPI_COMM_WORLD);
	
	char f_name[100];
	sprintf(f_name, "cliqueEdges_Processor_%d.mtx", my_rank);
	const char* fileName = f_name; // = "output.txt";  // Specify the file name		
	char filePath[256]; // Adjust the size as needed
	snprintf(filePath, sizeof(filePath), "%s/%s", folderName, fileName);
	FILE* tempFile = fopen(filePath, "a");
	

	// MPI_Barrier(MPI_COMM_WORLD);
	
	// Initilizing the adjMatrix
	if (my_rank == 0){
		printf("array size: %d \n", arraySize);
		int lp, mp , rp, edges;
		
		
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
		sscanf(line, "%d %d %d \n", &lp , &rp, &edges);
		printf("%d %d %d \n", lp, rp, edges);
		
		if (rp == 0){
			sscanf(line, "%d %d %d %d\n", &lp , &mp, &rp, &edges);
			printf("Found right partition to be zero.\n");
			printf("%d %d %d \n", lp, rp, edges);
		}
		
		for (int i = 0; i < edges; i++) {
			int row, col;
			fscanf(file, "%d %d \n", &row, &col); // , &temp
			adjMatrix_buffer[ ((col - adjOffset) * graphNodes) + (row - adjOffset) ] = 1;
			m_hat++;
		}
		InitialEdges = m_hat;
		printf("Initial edges: %d\n",m_hat);
		fclose(file);

		int* S  = (int*)malloc(graphNodes * sizeof(int));
		
		// char f_name[100];
		// sprintf(f_name, "cliqueEdgesProcessor_%d.mtx", my_rank);
		// const char* fileName = f_name; // = "output.txt";  // Specify the file name
		int* edgesRemoved = (int*)malloc(comm_size * sizeof(int));
		MPI_Request gRequest, sRequest;		
		int k_hat = get_k_hat(graphNodes, m_hat, delta);
		
		for ( int r = 1; r < comm_size; r++){
			MPI_Send(&k_hat, 1, MPI_INT, r, 0, MPI_COMM_WORLD);
		}
		int sendBufferSize = k_hat + 2;
		int* send_buf = (int*)malloc(sendBufferSize * sizeof(int));
		int q = graphNodes;
		srand(time(NULL));
		printf("k_hat: %d\n", k_hat);
		while (k_hat > 1) {
			for (int i  = 0; i < graphNodes; i++)
				S[i] = i;
			int setSize = graphNodes;
			int sendSize = k_hat + 2;
			while( setSize > 0){
				int R = comm_size;
				// printf("SetSize %d\n", setSize);
				int edgeRemoved = 0;
				int r = 0;
				int remainingVertices;

				for (r = 1; r < comm_size; r++){
					edgesRemoved[r] = 0;
					send_buf[1] = q;
					q++;
					if (setSize < k_hat)
						remainingVertices = setSize;
					else
						remainingVertices = k_hat;		
					send_buf[0] = remainingVertices;
					for (int k = 0; k < remainingVertices; k++){
						int randomVertex = rand() % setSize;
						send_buf[k + 2]	= S[randomVertex];			
						for (int j = randomVertex; j < setSize - 1; j++){
							S[j] = S[j + 1];
						}
						setSize--;						
					}	
					MPI_Send(send_buf, sendBufferSize, MPI_INT, r, 0, MPI_COMM_WORLD); // , &request			 , &sRequest	
					if ((setSize <= k_hat || r == comm_size - 1) && setSize > 0){
						R = r + 1;
						// printf("MPI process %d finds common neighbour for %d remaining vertices with %d processors.\n", my_rank, setSize, R -1);
						int remainingVertices;
						if (setSize < k_hat)
							remainingVertices = setSize;
						else
							remainingVertices = k_hat;
						
						int* targetIdx = (int*)malloc(remainingVertices * sizeof(int));
						int* neighbours = (int*)malloc(graphNodes * sizeof(int));
						int neighboursIdx = 0;
						printf("Processor %d with %d k_hat (", my_rank, remainingVertices);
						for(int k = 0; k < remainingVertices; k++){
							int randomVertex = rand() % setSize;
							// printf("%d ", S[randomVertex]);
							int p = (int) floor((S[randomVertex] * graphNodes)/arraySize);
							int idx = (S[randomVertex] * graphNodes) % arraySize;
							targetIdx[k] = (p - my_rank) * arraySize + idx;
							// printf("Idx %d\n", targetIdx[k]);
							for (int j = randomVertex; j < setSize - 1; j++){
								S[j] = S[j + 1];
							}
							setSize--;
						}
							
						printf(") \n");
						
						for(int i = 0; i < graphNodes; i++){
							int temp = 1;
							for(int k = 0; k < remainingVertices; k++){
								temp = temp * adjMatrix_buffer[targetIdx[k] + i];
							}
							if (temp){
								neighbours[neighboursIdx] = i;	
								edgeRemoved += remainingVertices;
								neighboursIdx++;
								for(int k = 0; k < remainingVertices; k++){
									adjMatrix_buffer[targetIdx[k] + i] = 0;
								}
							}					
						}
						for(int i = 0; i < neighboursIdx; i++){
							fprintf(tempFile, "%d %d\n", neighbours[i], q);
							edgesAddedToClique++;
						}
						for(int v = 0; v < remainingVertices; v++){
							fprintf(tempFile, "%d %d\n", q, S[v]);			
							edgesAddedToClique++;
						}
						// fclose(tempFile);
						free(targetIdx);
						free(neighbours);
						q++;
						printf("processor %d edges removed: %d\n", my_rank, edgeRemoved);
						break;
					}
					if (setSize == 0 ){
						R = r + 1;
						break;
					}
				}
				edgesRemoved[0] = edgeRemoved;
				m_hat = m_hat - edgesRemoved[0];
				// printf("p %d edges removed: %d\n", my_rank, edgesRemoved[0]);
				for ( int r = 1; r < R; r++){
					MPI_Recv(&edgeRemoved, 1, MPI_INT, r, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					edgesRemoved[r] = edgeRemoved;
					m_hat = m_hat - edgesRemoved[r];
					// printf("p %d edges removed: %d\n", r, edgesRemoved[r]);
				}
				// printf("remaining edges : %d\n", m_hat);
			}
			k_hat = get_k_hat(graphNodes, m_hat, delta);	
			printf("k_hat: %d\n", k_hat);
		}
		printf("Sending termination request\n");
		send_buf[0] = 1;
		for ( int r = 1; r < comm_size; r++){
			MPI_Send(send_buf, sendBufferSize, MPI_INT, r, 0, MPI_COMM_WORLD); // , &request			 , &sRequest		
		}
		
		// MPI_Wait(&sRequest, MPI_STATUS_IGNORE);
		TotalCliqueEdges = edgesAddedToClique;
		
		for ( int r = 1; r < comm_size; r++){
			MPI_Recv(&edgesAddedToClique, 1, MPI_INT, r, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			TotalCliqueEdges += edgesAddedToClique;
		}
		printf("Edges in the cliques: %d | Remaining/Trivial edges: %d | Edges in compressed graph: %d \n", TotalCliqueEdges, m_hat, TotalCliqueEdges + m_hat);
		free(send_buf);
		free(edgesRemoved);
		free(S);
	}
	else{
		MPI_Request request;
		int received;
		char f_name[100];
		

		int k_hat_buffer;
		MPI_Recv(&k_hat_buffer, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		
		int recvBufferSize = k_hat_buffer + 2;
		int* recv_buf = (int*)malloc(recvBufferSize * sizeof(int));
		while(1){
			
			MPI_Recv(recv_buf, recvBufferSize, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			int k_hat = recv_buf[0];
			int Q = recv_buf[1];
			if (k_hat < 2)
				break;
			printf("Receiving Processor %d waiting Line 255 with %d k_hat (", my_rank, k_hat);	
			for(int k = 0; k < k_hat; k++){
				printf("%d ", recv_buf[k + 2]);
			}
			printf(") \n");
			int edgeRemoved = 0;
			int* targetIdx = (int*)malloc(k_hat * sizeof(int));	
			int* neighbours = (int*)malloc((graphNodes) * sizeof(int));				
			for(int k = 0; k < k_hat; k++){
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
					neighboursIdx++;
					edgeRemoved += k_hat;
					for(int k = 0; k < k_hat; k++){
						adjMatrix_buffer[targetIdx[k] + i] = 0;
					}
				}
			}
			MPI_Send(&edgeRemoved, 1, MPI_INT, 0, 0, MPI_COMM_WORLD); // , &request
			printf("processor %d edges removed: %d\n", my_rank, edgeRemoved);
			
			for(int i = 0; i < neighboursIdx; i++){
				fprintf(tempFile, "%d %d\n", neighbours[i], Q);
				edgesAddedToClique++;
			}
			for(int v = 0; v < k_hat; v++){
				fprintf(tempFile, "%d %d\n", Q, recv_buf[v + 2]);			
				edgesAddedToClique++;
			}
			
			free(targetIdx);
			free(neighbours);
		}
		// printf("Processor %d terminating as (k_hat) < 2 \n", my_rank);
		MPI_Send(&edgesAddedToClique, 1, MPI_INT, 0, 0, MPI_COMM_WORLD); // , &request
		free(recv_buf);
	}

    
	
	// Writing Trivial edges to File
	for(int j = 0; j < arraySize; j++){
		if (adjMatrix_buffer[j] == 1){
			int realativeIdx = j + (my_rank * arraySize);
			int w = (int) floor (realativeIdx/graphNodes); // col index
			int u = (int) (realativeIdx % graphNodes);  // row index
			fprintf(tempFile, "%d %d\n", u, w);
		}	
	}
	
	MPI_Barrier(MPI_COMM_WORLD);

    // Destroy the window
    // printf("[MPI process %d] adjMatrix_window destroyed.\n", my_rank);
	MPI_Win_free(&adjMatrix_window);
	printf("[MPI process %d] adjMatrix_window destroyed.\n", my_rank);
    
	// MPI_Win_free(&whileTerminator);
	if (my_rank == 0){
		double end = MPI_Wtime();
		printf("Time elapsed during the job: %.2fs.\n", end - start);
		printf("Compression ratio: %f \n", (float) InitialEdges/(float) (TotalCliqueEdges + m_hat)); 
	}

	fclose(tempFile);
	printf("[MPI process %d] temp file destroyed.\n", my_rank);
	// MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize();
    return 0;
}

