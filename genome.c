#include "genome.h"

static int** load(char *filename, int *n, int *m, int** fs) {
	FILE *fp;
	if(!(fp = fopen(filename, "rb"))) {
		fprintf(stderr, "\nFile could not be opened.\n");
		return NULL;
	}

	// Read n
	if(fread(n, sizeof(int), 1, fp) == 0) {
		fprintf(stderr, "\nShort read error.\n");
	}

	// Read m
	if(fread(m, sizeof(int), 1, fp) == 0) {
		fprintf(stderr, "\nShort read error.\n");
	}

	int** mainArray = malloc(sizeof(*mainArray) * (*n));
	int* firstSeq = malloc(sizeof(*firstSeq) * (*n));

	int mTemp, nTemp, temp, count;

	// Read first sequence
	for(nTemp = 0; nTemp < *n; nTemp++) {
		if(fread(&(firstSeq[nTemp]), sizeof(int), 1, fp) == 0) {
			fprintf(stderr, "\nFailed to read first sequence.\n");
		}
	}

	// Move file pointer back
	fseek(fp, sizeof(int) * 2, SEEK_SET);

	for(nTemp = 0; nTemp < *n; nTemp++) {
		int* tempSeq = malloc(sizeof(*tempSeq) * (*m));
		mainArray[nTemp] = tempSeq;
	}

	// Load the 2D array with information
	for(mTemp = 0; mTemp < *m; mTemp++) {
		count = 0;
		for(nTemp = 0; nTemp < *n; nTemp++) {
			if(fread(&temp, sizeof(int), 1, fp) == 0) {
				fprintf(stderr, "\nShort read in sequences\n");
			}
			mainArray[temp - 1][mTemp] = count++;
		}
	}
	
	// Re-organize the arrays
	int** fixedArray = malloc(sizeof(*fixedArray) * (*n));
	for(nTemp = 0; nTemp < *n; nTemp++) {
		fixedArray[nTemp] = mainArray[firstSeq[nTemp] - 1];
	}
	free(mainArray);
	*fs = firstSeq;
	fclose(fp);
	return fixedArray;
}

typedef struct adjList {
	struct element *e;
	struct adjList *next;
} AdjList;

typedef struct element {
	int idx;
	int val;
	int *coord;
	struct adjList *list;
	int connections;
	int weight;
} Element;

static void connect(Element *a, Element *b) {
	AdjList *aList = a->list;
	a->list = malloc(sizeof(*(a->list)));
	AdjList newConn = {.e = b, .next = aList};
	*(a->list) = newConn;
	a->connections++;
}

static void initialize(Element *a, Element *b, int max) {
	int *a_coord = a->coord;
	int *b_coord = b->coord;
	int a_idx = 0, b_idx = 0;
	while(a_coord[a_idx] < b_coord[b_idx]) {
		if(a_idx == max - 1) {
			connect(a, b);
			break;
		}
		a_idx++;
		b_idx++;
	}
}

// Create adjacency list graph
static Element* makeGraph(int* firstSeq, int** coordinates, int n, int m) {
	int nTemp, nTemp2;

	// Initialize elements
	Element* elementArray = malloc(sizeof(*elementArray) * n);
	for(nTemp = 0; nTemp < n; nTemp++) {
		Element newElement = {.idx = nTemp, .val = firstSeq[nTemp], .coord = coordinates[nTemp], .list = NULL, .connections = 0, .weight = 0};
		elementArray[nTemp] = newElement;
	}

	// Free coordinate array and first sequence
	free(coordinates);
	free(firstSeq);

	for(nTemp = 0; nTemp < n - 1; nTemp++) {
		for(nTemp2 = nTemp; nTemp2 < n; nTemp2++) {
			initialize(&(elementArray[nTemp]), &(elementArray[nTemp2]), m);
		}
	}
	return elementArray;
}

static void getSequence(int* seq, Element *eList, int seqIdx, int idx) {
	seq[seqIdx] = eList[idx].val;
	if(eList[idx].weight == 1) {
		return;
	}

	AdjList *aList = eList[idx].list;
	while(aList != NULL) {
		if(aList->e->weight != eList[idx].weight - 1) {
			aList = aList->next;
		} else {
			getSequence(seq, eList, seqIdx + 1, aList->e->idx);
			break;
		}
	}

}
static void auxDFS(Element *eList, int idx) {
	AdjList *aList = eList[idx].list;
	if(aList == NULL) {
		eList[idx].weight = 1;
		return;
	}
	int largest = 0;
	while(aList != NULL) {
		if(aList->e->weight > largest) {
			largest = aList->e->weight;
		}
		aList = aList->next;
	}
	eList[idx].weight = largest + 1;
}

static void getLengthDFS(Element *eList, int idx) {
	int i;
	for(i = idx; i >= 0; i--) {
		if(eList[i].weight == 0) {
			auxDFS(eList, i);
		}
	}
}
static int* longestSequence(Element *eList, int *size, int n) {
	int nTemp, longest = 0, startIdx = 0;

	//Find the start index of longest sequence
	for(nTemp = 0; nTemp < n; nTemp++) {
		if(eList[nTemp].connections > longest) {
			startIdx = nTemp;
			longest = eList[nTemp].connections;
		}
	}

	int store = 0;
	getLengthDFS(eList, n -1);
	store = eList[startIdx].weight;
	*size = store;

	int* seq = malloc(sizeof(*seq) * store);

	getSequence(seq, eList, 0, startIdx);
	return seq;
}

static void freeElements(Element *eList, int n) {
	int i;
	for(i = 0; i < n; i++) {
		free(eList[i].coord);
		while(eList[i].list != NULL) {
			AdjList *temp = eList[i].list;
			eList[i].list = eList[i].list->next;
			free(temp);
		}
	}
	free(eList);
}

int *Longest_conserved_gene_sequence(char *filename, int *size_of_seq) {
	int n = 0, m = 0;
	int* firstSeq = NULL;

	// Generate multi dimensional array in 2D coordinate represenation
	int** sequenceMap = load(filename, &n, &m, &firstSeq);
	if(sequenceMap == NULL) {
		*size_of_seq = 0;
		return NULL;
	}

	// Generate graph using element struct type
 	Element *eList = makeGraph(firstSeq, sequenceMap, n, m);

	// Traverse graph to find longest sequence and its length
	int size = 0;
	int* ans = longestSequence(eList, &size, n);
	*size_of_seq = size;

	// Clean memory
	freeElements(eList, n);
	return ans;
}

