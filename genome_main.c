#include "genome.h"

int main(int argc, char **argv) {
	char* inputFileName = NULL;

	if(argc != 2) {
		fprintf(stderr, "Non-standard number of arguments.");
	} else {
		inputFileName = argv[1];
	}
	int seqSize = 0;
	int *seq = Longest_conserved_gene_sequence(inputFileName, &seqSize);

	free(seq);
	fprintf(stdout, "Length: %d\n", seqSize);
	return EXIT_SUCCESS;
}
