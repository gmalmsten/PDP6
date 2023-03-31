#include "stencil.h"
#include <string.h>


int main(int argc, char **argv) {
	if (4 != argc) {
		printf("Usage: stencil input_file output_file number_of_applications\n");
		return 1;
	}
	char *input_name = argv[1];
	char *output_name = argv[2];
	int num_steps = atoi(argv[3]);

	// Read input file
	double *input;
	int num_values;
	if (0 > (num_values = read_input(input_name, &input))) {
		return 2;
	}

	// Stencil values
	double h = 2.0*PI/num_values;
	const int STENCIL_WIDTH = 5;
	const int EXTENT = STENCIL_WIDTH/2;
	const double STENCIL[] = {1.0/(12*h), -8.0/(12*h), 0.0, 8.0/(12*h), -1.0/(12*h)};

	double input_test[20];
	for(int i = 0; i< 20; i++) input_test[i] = (double)i;

	// Start timer
	double start = MPI_Wtime();
	
	// Initialize MPI and assign sub lists
	MPI_Init(&argc, &argv);
	int rank, right, left, num_proc;
	MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int chunkSz = num_values/num_proc;
	double * sub_list = (double *)malloc(chunkSz*sizeof(double));
	memcpy(sub_list, &input[rank*chunkSz], chunkSz*sizeof(double));

	// Create circular topology
	MPI_Comm CIRC_COMM;
	int dims[1];
	int periods[1];
	dims[0] = num_proc;
	periods[0] = 1;
	int reorder = 0;
	MPI_Cart_create(MPI_COMM_WORLD, 1, dims, periods, reorder, &CIRC_COMM);
	MPI_Comm_rank(CIRC_COMM, &rank);
	MPI_Cart_shift(CIRC_COMM, 0, -1, &right, &left);


	// Allocate data for result
	double *output;
	if (NULL == (output = malloc(num_values * sizeof(double)))) {
		perror("Couldn't allocate memory for output");
		return 2;
	}
	// Repeatedly apply stencil
	for (int s=0; s<num_steps; s++) {
		// Apply stencil
		for (int i=0; i<EXTENT; i++) {
			double result = 0;
			for (int j=0; j<STENCIL_WIDTH; j++) {
				int index = (i - EXTENT + j + num_values) % num_values;
				result += STENCIL[j] * input[index];
			}
			output[i] = result;
		}
		for (int i=EXTENT; i<num_values-EXTENT; i++) {
			double result = 0;
			for (int j=0; j<STENCIL_WIDTH; j++) {
				int index = i - EXTENT + j;
				result += STENCIL[j] * input[index];
			}
			output[i] = result;
		}
		for (int i=num_values-EXTENT; i<num_values; i++) {
			double result = 0;
			for (int j=0; j<STENCIL_WIDTH; j++) {
				int index = (i - EXTENT + j) % num_values;
				result += STENCIL[j] * input[index];
			}
			output[i] = result;
		}
		// Swap input and output
		if (s < num_steps-1) {
			double *tmp = input;
			input = output;
			output = tmp;
		}
	}
	free(input);
	MPI_Finalize();
	// Stop timer
	double my_execution_time = MPI_Wtime() - start;

	// Write result
	printf("Took %fs\n", my_execution_time);
#ifdef PRODUCE_OUTPUT_FILE
	if (0 != write_output(output_name, output, num_values)) {
		return 2;
	}
#endif

	// Clean up
	free(output);

	return 0;
}


int read_input(const char *file_name, double **values) {
	FILE *file;
	if (NULL == (file = fopen(file_name, "r"))) {
		perror("Couldn't open input file");
		return -1;
	}
	int num_values;
	if (EOF == fscanf(file, "%d", &num_values)) {
		perror("Couldn't read element count from input file");
		return -1;
	}
	if (NULL == (*values = malloc(num_values * sizeof(double)))) {
		perror("Couldn't allocate memory for input");
		return -1;
	}
	for (int i=0; i<num_values; i++) {
		if (EOF == fscanf(file, "%lf", &((*values)[i]))) {
			perror("Couldn't read elements from input file");
			return -1;
		}
	}
	if (0 != fclose(file)){
		perror("Warning: couldn't close input file");
	}
	return num_values;
}


int write_output(char *file_name, const double *output, int num_values) {
	FILE *file;
	if (NULL == (file = fopen(file_name, "w"))) {
		perror("Couldn't open output file");
		return -1;
	}
	for (int i = 0; i < num_values-1; i++) {
		if (0 > fprintf(file, "%.4f ", output[i])) {
			perror("Couldn't write to output file");
		}
	}
	if (0 > fprintf(file, "%.4f", output[num_values])){
		perror("Couldn't write to output file");
	}
	if (0 > fprintf(file, "\n")) {
		perror("Couldn't write to output file");
	}
	if (0 != fclose(file)) {
		perror("Warning: couldn't close output file");
	}
	return 0;
}
