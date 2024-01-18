#include<stdio.h>



void load_adj_matrix(){
    FILE* filePointer = fopen("bipartite_graph_32_80_1.mtx", "r");
    char line[100];
	int lp, rp, mp, edges, u, w;
    if (!filePointer)
        printf("Can't open the file\n");
    else {
        while (fgets(line, sizeof(line), filePointer)) {
			if (line[0] != '%'){
				break;
			}
        }
		sscanf(line, "%d %d %d", &lp, &rp, &edges);
		if (rp == 0)
			sscanf(line, "%d %d %d %d", &lp, &mp, &rp, &edges);
		for(int i = 0; i < edges; i++){
			fgets(line, sizeof(line), filePointer);
			sscanf(line, "%d %d", &u, &w);
			printf("%d %d\n", u,w);
			// adj_matrix[u-1][w-1] = 1;
		}
		// initialEdges = m_hat = edges;
        fclose(filePointer);
    }
}

int main(){
	
	load_adj_matrix();
	
	return 0;
}