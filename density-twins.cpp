/* This program computes all pairs of density twins of order V for all V=2..9
* By the theorem in the paper "density twins", we know that density twins must have the same order.
* Furthermore, there are no density twins of order less than 7.
* Running the program for V=7 yields a unique pair of density twins and for V=8 yields four pairs of density twins.
*/

#include <stdio.h>

/* To run the program, it is necessary to use a formatted file containing the adjacency matrix of all graphs on V vertices, with the name graphs<V>.txt (e.g., graphs8.txt).
This file (for any particaular V=2..9) is available from the author, or can also be used from McKay's database:  https://users.cecs.anu.edu.au/~bdm/data/graphs.html after
expanding it from the dense format specified there, using the instructions given on that page. */
#define DATAFILE "d:\\research\\general\\combinatorial data\\graphs"

const int V = 8; /* plug here any value in the range 1..9 to run the program; V=8 runs in less than an hour, V=9 is lengthy. */

const int numGraphs[] = { 1, 1, 2, 4, 11, 34, 156, 1044, 12346, 274668 };  //number of graphs on V vertices for V=0,...,9.
short graph[274668][V][V];  //graph[r] rpresents the adjacency matrix of the r'th graph on V vertices
short homog[274668][4]; // for graphs[r], homog[r][0] = independence number, homog[r][1] = clique number, homog [r][2] = #of empty twins, homog[r][3] = # of connected twins.
int homomorph[V];  // represents a homomorphism from some graph to some other graph.(i.e., a function from 0..V-1 to 0..V-1)
int loopVec[V]; // represents the binary vector for the possible vertices with loops around them.



// load all graphs of order n to graph
void loadGraph(int n) {
	FILE* datafile = NULL;
	char dataFileName[200];
	char buffer[100];
	sprintf_s(dataFileName, "%s%d.txt", DATAFILE, n);
	fopen_s(&datafile, dataFileName, "r");
	for (int r = 0; r < numGraphs[n]; r++) {
		fgets(buffer, sizeof(buffer), datafile);
		fgets(buffer, sizeof(buffer), datafile);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++)
				graph[r][i][j] = fgetc(datafile) - '0';
			fgetc(datafile);
		}
	}
	fclose(datafile);
}

// print the adjacency matrix of graph[r] on n vertices
void printGraph(int n, int r) {
	printf("Printing graph #%d on %d vertices\n", r, n);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++)
			printf("%d", graph[r][i][j]);
		printf("\n");
	}
}

// count the number of empty twins and connected twins in graph[r] on n vertices
void countTwins(int n, int r) {
	int twinConn = 0;
	int twin = 0;
	for (int i = 0; i < n; i++)
		for (int j = i + 1; j < n; j++) {
			bool found = true;
			for (int k = 0; k < n; k++)
				if (graph[r][i][k] != graph[r][j][k] && k != i && k != j) {
					found = false;
					break;
				}
			if (!found)
				continue;
			graph[r][i][j] == 1 ? twinConn++ : twin++;
		}
	homog[r][2] = twin;
	homog[r][3] = twinConn;
}

// check if the current subset of vertices designated by loopVec is a clique (if edge == 1) or independent set (if edge == 0)
bool checkHomog(int n, int r, int edge) {
	for (int i = 0; i < n; i++)
		for (int j = i + 1; j < n; j++)
			if (loopVec[i] == 1 && loopVec[j] == 1 && graph[r][i][j] != edge)
				return false;
	return true;
}

// if edge==0, compute the independence number, if edge==1 compute the clique number
void maxHomog(int n, int r, int edge) {
	homog[r][edge] = 0;
	bool first = true;
	for (int i = 0; i < n; i++)
		loopVec[i] = 1;
	while (true) {
		int j = 0;
		while (loopVec[j] == 1 && j < n) {
			loopVec[j] = 0;
			j++;
		}
		if (j == n && first)
			first = false;
		else if (j == n)
			return;
		else
			loopVec[j] = 1;
		int size = 0;
		for (int i = 0; i < n; i++)
			if (loopVec[i] == 1)
				size++;
		if (size > homog[r][edge] && checkHomog(n, r, edge))
			homog[r][edge] = size;
	}
}

// check if there is a current loop homomorphism from graph[r1] on n1 vertices to graph[r2], both on n vertices
bool checkLoopVecHom(int n, int r1, int r2) {
	for (int i = 1; i < n; i++)
		for (int j = i + 1; j < n; j++) {
			if (homomorph[i] != homomorph[j] && graph[r1][i][j] != graph[r2][homomorph[i]][homomorph[j]])
				return false;
			if (homomorph[i] == homomorph[j] && (graph[r1][i][j] != loopVec[homomorph[i]]))
				return false;
		}
	return true;
}

// just like the previous function but with a ``speed hack''
bool checkLoopVecHomRand(int n, int r1, int r2) {
	for (int i = 0; i < n; i += 2)
		for (int j = i + 1; j < n; j++) {
			if (homomorph[i] != homomorph[j] && graph[r1][i][j] != graph[r2][homomorph[i]][homomorph[j]])
				return false;
			if (homomorph[i] == homomorph[j] && (graph[r1][i][j] != loopVec[homomorph[i]]))
				return false;
		}
	return true;
}

// check if there is any homomorphism for the current loopVec
bool checkLoopVec(int n, int r1, int r2) {
	for (int i = 0; i < n; i++)
		homomorph[i] = 0;
	if (checkLoopVecHom(n, r1, r2))  // check trivial homomorphism and current loopvec
		return true;
	while (true) {
		int j = 0;
		while (homomorph[j] == n - 1 && j < n) {
			homomorph[j] = 0;
			j++;
		}
		if (j == n)
			return false; // all homomorphisms for the current loopVec are invalid, so return false
		homomorph[j]++;
		if (checkLoopVecHomRand(n, r1, r2) && checkLoopVecHom(n, r1, r2))
			return true;
	}
}

bool checkDeep(int n, int r1, int r2, int skip) {
	int counter = 0;
	bool first = true;
	for (int i = 0; i < n; i++)
		loopVec[i] = 1;
	while (true) {
		counter++;
		int j = 0;
		while (loopVec[j] == 1 && j < n) {
			loopVec[j] = 0;
			j++;
		}
		if (j == n && first)
			first = false;
		else if (j == n)
			return true;
		else
			loopVec[j] = 1;
		if (counter % skip != 0)
			continue;
		if (!checkLoopVec(n, r1, r2))
			return false;
	}
}

// check if grpah[s1][r1] and graph[s2[r2] are compatible
bool checkCompatible(int n, int r1, int r2) {
	for (int i = 0; i < n; i++)
		loopVec[i] = 0;
	if (!checkLoopVec(n, r1, r2) || !checkLoopVec(n, r2, r1))
		return false;
	for (int i = 0; i < n; i++)
		loopVec[i] = 1;
	if (!checkLoopVec(n, r1, r2) || !checkLoopVec(n, r2, r1))
		return false;
	if (!checkDeep(n, r1, r2, 81) || !checkDeep(n, r2, r1, 79)) {
		return false;
	}
	if (!checkDeep(n, r1, r2, 27) || !checkDeep(n, r2, r1, 25)) {
		return false;
	}
	if (!checkDeep(n, r1, r2, 9) || !checkDeep(n, r2, r1, 7)) {
		return false;
	}
	if (!checkDeep(n, r1, r2, 3) || !checkDeep(n, r2, r1, 3)) {
		return false;
	}
	// Reached brute force
	return checkDeep(n, r1, r2, 1) && checkDeep(n, r2, r1, 1);
}

int main()
{
	loadGraph(V);

	for (int r = 0; r < numGraphs[V]; r++) {
		maxHomog(V, r, 0);
		maxHomog(V, r, 1);
		countTwins(V, r);
	}

	for (int r1 = 0; r1 < numGraphs[V]; r1++) {
		if (r1 % 100 == 0)
			printf("%d\n", r1);
		if (homog[r1][2] == 0 || homog[r1][3] == 0) // graph[r1] is denisty atomic if it has no full twins or no empty twins
			continue;
		for (int r2 = r1 + 1; r2 < numGraphs[V]; r2++) {
			if (homog[r1][2] != homog[r2][2] || homog[r1][3] != homog[r2][3] || homog[r1][0] != homog[r2][0] || homog[r1][1] != homog[r2][1])
				continue;  // graph[r1] and graph[r2] cannot be compatible
			if (checkCompatible(V, r1, r2)) {
				printf("Graph %d and graph %d on %d vertices are compatible\n", r1, r2, V);
				printf("\n");
				printGraph(V, r1);
				printf("\n");
				printGraph(V, r2);
			}
		}
	}
}