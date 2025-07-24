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
    if (de < 10){
        de = 10;
    }
    float nu = delta * log2((double)graph_nodes);
    k_hat = floor((double)nu / (log2((double)de)));
	return k_hat;
}

int main(int argc, char* argv[]) {
	const char* folderName = "temp"; // Specify the folder name
	const char* compressedGraphFolderName = "restrecturedGraph"; // Specify the folder name
	createFolder(folderName);
	createFolder(compressedGraphFolderName);
    MPI_Init(&argc, &argv);
	double start = MPI_Wtime();
    double saveCliqueStart, saveCliqueTime, mergeFilesStart, mergeFilesTime, readingTimeStart, readingTimeEnd;
    double totalWriteTime = 0.0;
	saveCliqueTime = mergeFilesTime = 0.0;
    int comm_size;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    float weight = 1.0; // default weight for unweighted graphs
	const char* inputFilename = argv[1];
	
	float delta = atof(argv[2]);
	int density = 99; //atoi(argv[4]);
	int instance = atoi(argv[3]);
    
  
    // Get my rank
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	
	int edgesAddedToClique = 0;
	int TotalCliqueEdges = 0;
	int InitialEdges = 0;
	int m_hat = 0;
	int total_cliques;
    // Create the window for a 10x10 matrix

    bool* adjMatrix_buffer;
	
	
	int lp = 0, rp = 0, edges = 0;

    // Rank 0 reads the file
    if (my_rank == 0) {
        FILE* inputFile = fopen(inputFilename, "r");
        if (!inputFile) {
            perror("Error opening file");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        char line[1024];
        while (fgets(line, sizeof(line), inputFile)) {
            if (line[0] == '%') continue;
            if (sscanf(line, "%d %d %d", &lp, &rp, &edges) == 3) break;
        }

        fclose(inputFile);
    }

    // Broadcast lp, rp, and edges to all ranks
    MPI_Bcast(&lp, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&rp, 1, MPI_INT, 0, MPI_COMM_WORLD);
    // MPI_Bcast(&edges, 1, MPI_INT, 0, MPI_COMM_WORLD);

	int graphNodes = (lp > rp) ? lp : rp; // = atoi(argv[2])
	int arraySize = (int) ceil((double)(graphNodes * graphNodes)/(double)comm_size);
    const int ROWS = graphNodes;
    const int COLS = graphNodes/comm_size;

    MPI_Win adjMatrix_window;
    MPI_Win_allocate_shared(arraySize * sizeof(bool), sizeof(bool), MPI_INFO_NULL, MPI_COMM_WORLD, &adjMatrix_buffer, &adjMatrix_window);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_File readFile;
    MPI_Offset filesize, chunksize, offset;
	// int offset;
    MPI_File_open(MPI_COMM_WORLD, inputFilename, MPI_MODE_RDONLY, MPI_INFO_NULL, &readFile);
    MPI_File_get_size(readFile, &filesize);
    chunksize = filesize / comm_size;
	char *buffer; // store the block of file read by each processor.
	int adjOffset = 1;
	char f_name[100];
	sprintf(f_name, "p_%d.mtx", my_rank);
	const char* tempfileName = f_name; // = "output.txt";  // Specify the file name		
	char filePath[256]; // Adjust the size as needed
	snprintf(filePath, sizeof(filePath), "%s/%s", folderName, tempfileName);
	FILE* tempFile = fopen(filePath, "a");
	MPI_Request request;
	// // // Reading the file in parallel
	
	if (my_rank == 0){
		//printf("%s, %d, %f, %d, %d\n", inputFilename, graphNodes, delta, density, instance);
		readingTimeStart = MPI_Wtime();
		int lp, mp , rp, edges;
		offset = 0;
		buffer = (char *)malloc((chunksize) * sizeof(char));
		MPI_File_read_at(readFile, offset, buffer, chunksize, MPI_BYTE, MPI_STATUS_IGNORE);
		int i = 0;
		int lastChar = buffer[chunksize - i];
		while (lastChar !=  10 ) {
			int lastChar = buffer[chunksize - i];
			// printf("Character %d at idx %ld\n", lastChar, chunksize - i);
			if (lastChar ==  10){
				break;
			}
			i++;
		}
		offset = offset + chunksize - i;
		//printf("offset is %ld\n", offset );
		MPI_Send(&offset, 1,MPI_LONG_LONG , my_rank + 1, 0, MPI_COMM_WORLD);
		char *line = strtok(buffer, "\n");
		int length = line ? line - buffer : 0;
		while (line[0] ==  '%'){
			line = strtok(NULL, "\n");
		}
		sscanf(line, "%d %d %d \n", &lp , &rp, &edges);
		
		if (rp == 0){
			sscanf(line, "%d %d %d %d\n", &lp , &mp, &rp, &edges);
		}
		// printf("Process %d: %d %d %d\n", my_rank, lp, rp, edges);
		InitialEdges = edges;
		m_hat = edges;
		line = strtok(NULL, "\n");
		while (line != NULL && length <= chunksize - i) {
			// Process the line as needed
			int row, col;
			sscanf(line, "%d %d \n", &row, &col); // , &temp
			adjMatrix_buffer[ ((col - adjOffset) * graphNodes) + (row - adjOffset) ] = 1;
			// printf("Process %d (%d): %d %d %d\n", my_rank, ((col - adjOffset) * graphNodes) + (row - adjOffset) , row, col, adjMatrix_buffer[((col - adjOffset) * graphNodes) + (row - adjOffset) ]);
			line = strtok(NULL, "\n");
			length = line ? line - buffer : 0;
			// printf("Length: %d \n", length);
		}
		
	}
	else{
		MPI_Recv(&offset, 1, MPI_LONG_LONG, my_rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		//printf("[%d] Offset %d \n", my_rank, offset);
		// chunksize += (my_rank * chunksize) - offset;
		if (my_rank == comm_size - 1) 
			chunksize = filesize - offset;
		buffer = (char *)malloc((chunksize) * sizeof(char));
		MPI_File_read_at(readFile, offset, buffer, chunksize, MPI_BYTE, MPI_STATUS_IGNORE);
		int i = 0;
		int lastChar = buffer[chunksize - i];

		while (lastChar !=  10 ){
			int lastChar = buffer[chunksize - i];
			// printf("Character %d at idx %ld\n", lastChar, chunksize - i);
			if (lastChar ==  10){
				break;
			}
			i++;
		}
		offset = offset + chunksize - i;
		// printf("chunksize is %ld\n", chunksize - i);
		if (my_rank != comm_size - 1){
			MPI_Send(&offset, 1, MPI_LONG_LONG, my_rank + 1, 0, MPI_COMM_WORLD);
			// printf("Processor [%d] sending to %d", my_rank, my_rank+1);
		}
		char *line = strtok(buffer, "\n");
		int length = line ? line - buffer : 0;
		while (line != NULL && length <= chunksize - i) {
			int row, col;
			sscanf(line, "%d %d \n", &row, &col); // , &temp
			int pos = ((col - adjOffset) * graphNodes) + (row - adjOffset);
			int p = (int) floor(pos/arraySize);
			int idx = pos % arraySize;
			int relIdx = (p - my_rank) * arraySize + idx;
			adjMatrix_buffer[relIdx] = 1;
			// printf("Process %d (%d): %d %d %d\n", my_rank, relIdx, row, col, adjMatrix_buffer[relIdx]);
			line = strtok(NULL, "\n");
			length = line ? line - buffer : 0;
			// printf("Length: %d \n", length);
		}
		
	}
	


	MPI_Barrier(MPI_COMM_WORLD);
	if (my_rank == 0){
		readingTimeEnd = MPI_Wtime() - readingTimeStart;
		//printf("Processor [%d] file read time: %lf s\n\n\n", my_rank, readingTimeEnd);
	}
	// if(my_rank == 0){
		// for(int i = 1; i <= graphNodes*graphNodes; i++){
			// printf("(%d) %d ",i-1, adjMatrix_buffer[i-1]);
			// if (i % graphNodes == 0)
				// printf("\n");
		// }
	// }
	// printf("\n");
	
	MPI_File_close(&readFile);
	MPI_Barrier(MPI_COMM_WORLD);
	// Initilizing the adjMatrix
	if (my_rank == 0){
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
		double* writeTimeBuffer = (double*)malloc(comm_size* sizeof(double));
		

		srand(time(NULL));
		//printf("k_hat: %d\n", k_hat);
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
						//printf("Processor %d with %d k_hat (", my_rank, remainingVertices);
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
							
						//printf(") \n");
						
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
						saveCliqueStart = MPI_Wtime();
						for(int i = 0; i < neighboursIdx; i++){
							fprintf(tempFile, "%d %d\n", neighbours[i], q);
							edgesAddedToClique++;
						}
						for(int v = 0; v < remainingVertices; v++){
							fprintf(tempFile, "%d %d\n", q, S[v]);			
							edgesAddedToClique++;
						}
						// fclose(tempFile);
						saveCliqueTime += (MPI_Wtime() - saveCliqueStart);
						free(targetIdx);
						free(neighbours);
						q++;
						//printf("processor %d edges removed: %d\n", my_rank, edgeRemoved);
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
			//printf("k_hat: %d\n", k_hat);
		}
		//printf("Sending termination request\n");
		printf("fnal Q:%d rank:%d\n", q - 1, my_rank);
		total_cliques = q-1;
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
		//printf("Edges in the cliques: %d | Remaining/Trivial edges: %d | Edges in restructured graph: %d \n", TotalCliqueEdges, m_hat, TotalCliqueEdges + m_hat);
		free(send_buf);
		free(edgesRemoved);
		free(S);
		
		saveCliqueStart = MPI_Wtime();
		for(int j = 0; j < arraySize; j++){
			if (adjMatrix_buffer[j] == 1){
				int realativeIdx = j + (my_rank * arraySize);
				int w = (int) floor (realativeIdx/graphNodes); // col index
				int u = (int) (realativeIdx % graphNodes);  // row index
				fprintf(tempFile, "%d %d\n", u, w);
			}	
		}
		saveCliqueTime += (MPI_Wtime() - saveCliqueStart);
		MPI_Gather(&saveCliqueTime, 1, MPI_DOUBLE, writeTimeBuffer , 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		
		for(int p = 0; p < comm_size; p++){
			totalWriteTime += writeTimeBuffer[p];
		}
		//printf("Total write Time = %lf", totalWriteTime );
		free(writeTimeBuffer);

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
			// printf("Q:%d rank:%d\n", Q, my_rank);
			if (k_hat < 2)
				break;
			//printf("Receiving Processor %d waiting Line 255 with %d k_hat (", my_rank, k_hat);	
			//for(int k = 0; k < k_hat; k++){
			//	printf("%d ", recv_buf[k + 2]);
			//}
			//printf(") \n");
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
			//printf("processor %d edges removed: %d\n", my_rank, edgeRemoved);
			saveCliqueStart = MPI_Wtime();
			for(int i = 0; i < neighboursIdx; i++){
				fprintf(tempFile, "%d %d\n", neighbours[i], Q);
				edgesAddedToClique++;
			}
			for(int v = 0; v < k_hat; v++){
				fprintf(tempFile, "%d %d\n", Q, recv_buf[v + 2]);			
				edgesAddedToClique++;
			}
			saveCliqueTime += (MPI_Wtime() - saveCliqueStart);
			free(targetIdx);
			free(neighbours);
		}
		// printf("Processor %d terminating as (k_hat) < 2 \n", my_rank);
		MPI_Send(&edgesAddedToClique, 1, MPI_INT, 0, 0, MPI_COMM_WORLD); // , &request
		
		saveCliqueStart = MPI_Wtime();
		for(int j = 0; j < arraySize; j++){
			if (adjMatrix_buffer[j] == 1){
				int realativeIdx = j + (my_rank * arraySize);
				int w = (int) floor (realativeIdx/graphNodes); // col index
				int u = (int) (realativeIdx % graphNodes);  // row index
				fprintf(tempFile, "%d %d\n", u, w);
			}	
		}
		saveCliqueTime += (MPI_Wtime() - saveCliqueStart);
		MPI_Gather(&saveCliqueTime, 1, MPI_DOUBLE, NULL, 0, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		free(recv_buf);
		

	}
	
	
	
	MPI_Barrier(MPI_COMM_WORLD);

    // Destroy the window
    // printf("[MPI process %d] adjMatrix_window destroyed.\n", my_rank);
	MPI_Win_free(&adjMatrix_window);
	//printf("[MPI process %d] adjMatrix_window destroyed.\n", my_rank);
    
	// MPI_Win_free(&whileTerminator);
	// if (my_rank == 0){
		// double end = MPI_Wtime();
		// printf("Time elapsed during the job: %.2fs.\n", end - start);
		// printf("Edge Reduction ratio: %f \n", (float) InitialEdges/(float) (TotalCliqueEdges + m_hat));
		// printf("%d, %.1f, %f, %.4f, %d\n", graphNodes, delta, (float) InitialEdges/(float) (TotalCliqueEdges + m_hat), end - start, comm_size);
	// }

	fclose(tempFile);
	//printf("[MPI process %d] temp file destroyed.\n", my_rank);
	// MPI_Barrier(MPI_COMM_WORLD);
	
	
	if (my_rank == 0){
		// printf("total cliques:%d rank:%d\n", total_cliques-graphNodes, my_rank);
		mergeFilesStart = MPI_Wtime();
		char f_name[100];
		char compressedFileName[256];
		char command[2400] = "cat "; // Adjust the size as needed
		sprintf(compressedFileName, "%s/restructured_graph_%d_%d_%d_%.1f.mtx", compressedGraphFolderName, graphNodes, density, instance, delta);
		FILE* outputFile = fopen(compressedFileName, "w");
		fprintf(outputFile, "%%MatrixMarket matrix coordinate pattern general\n");
		fprintf(outputFile, "%% Resulted restructured graph for given bipartite graph with %d nodes, %d density and %.1f delta.\n", graphNodes, density, delta);
		fprintf(outputFile, "%d %d %d %d\n", graphNodes,  total_cliques-graphNodes, graphNodes, TotalCliqueEdges+m_hat);
		fclose(outputFile);
		// strcat(command, "cat ");
		
		// snprintf(command, sizeof(command), "cat ", f_name, compressedFileName);
		for(int p = 0; p < comm_size; p++){
			sprintf(f_name, "%s/p_%d.mtx", folderName,p);
			strcat(command, f_name);
			strcat(command, " ");
			// snprintf(command, sizeof(command), "cat %s >> %s", f_name, compressedFileName);
			// int returnCode = system(command);
			// remove(f_name);
		}
		strcat(command, " >> ");
		strcat(command, compressedFileName);
		// printf("%s \n", command);
		int returnCode = system(command);
		mergeFilesTime += (MPI_Wtime() - mergeFilesStart);
		double end = MPI_Wtime();
		for(int p = 0; p < comm_size; p++){
			sprintf(f_name, "%s/p_%d.mtx", folderName,p);
			remove(f_name);
		}
		//system("rm -r temp");
		//printf("Processor [%d] file merging exe time: %lf s\n", my_rank, mergeFilesTime);
		printf("%d, %d, %d, %.1f, %f, %.4f, %.4f, %.4f, %.4f, %d\n", graphNodes, density, instance, delta, (float) InitialEdges/(float) (TotalCliqueEdges + m_hat), end - start, readingTimeEnd, totalWriteTime , mergeFilesTime, comm_size);
	}
	// printf("Processor [%d] file write time: %lf s\n", my_rank, saveCliqueTime);
	// double readingTimeStart = MPI_Wtime();
	MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    return 0;
}

