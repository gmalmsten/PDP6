#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

int read_input(const char *file_name, double **input){
    	FILE *file;
	if (NULL == (file = fopen(file_name, "r"))) {
		perror("Couldn't open input file");
		return -1;
	}

	int n;
    if (EOF == fscanf(file, "%d", &n)) {
        perror("Couldn't read element count from input file");
        return -1;
    }

    if (NULL == (*input = malloc(2 * n * n * sizeof(double)))) {
		perror("Couldn't allocate memory for matrix A");
		return -1;
	}


    for (int i=0; i<2*n*n; i++) {
		if (EOF == fscanf(file, "%lf", &((*input)[i]))) {
			perror("Couldn't read elements from input file");
			return -1;
		}
	}


    if (0 != fclose(file)){
		perror("Warning: couldn't close input file");
	}
	return n;
}

int main(int argc, char *argv[]){

    if (argc != 3) {
		printf("Usage: matmul input_file output_file\n");
		return 1;
	}
	char *input_name = argv[1];
	char *output_name = argv[2];

    // Initialize MPI
	MPI_Init(&argc, &argv);

	int rank, num_proc;
	MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Create cartesian grid virtual topology
	MPI_Comm GRID_COMM;
	int dims[2];
	int periods[2];
	int reorder = 0;  // May use reorder and MPI_Cart_rank(GRID_COMM, coords[], &rank);
	dims[0] = sqrt(num_proc);
    dims[1] = sqrt(num_proc);
	periods[0] = 0;
    periods[1] = 0;

	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &GRID_COMM);

    // Check that the processes have the correct coords
    int coords[2];
    MPI_Cart_coords(GRID_COMM, rank, 2, coords);
    printf("Rank %d coords (%d, %d)\n", rank, coords[0], coords[1]);
    // Initialize matrices
    double *A, *input;
    int n;

    // Source process reads input
    if(rank == 0){
    n = read_input(input_name, &input);
        printf("A: \n");

    // Print read input
    for(int i = 0; i<n; i++){
        for(int j = 0; j<n; j++){
            printf("%lf ", input[n*i + j]);
        }
        printf("\n");
    }

    printf("\n");
    printf("B: \n");
    for(int i = 0; i<n; i++){
        for(int j = 0; j<n; j++){
            printf("%lf ", input[n*n + n*i + j]);
        }
        printf("\n");
    }
    }

    // Broadcast set up parameters
    MPI_Bcast(&n, 1, MPI_INT, 0, GRID_COMM);


    // Define matrix type
    MPI_Datatype matrixtype;
    MPI_Type_vector(n/dims[0], n/dims[0], n, MPI_DOUBLE, &matrixtype);
    MPI_Type_commit(&matrixtype);


    // Scatter matrices to processes
    int sub_n = n/sqrt(num_proc);
    A = (double *)malloc(sub_n*sub_n*sizeof(double));
    MPI_Scatter(input, 1, matrixtype, A, sub_n*sub_n, MPI_DOUBLE, 0, GRID_COMM);



    printf("Rank %d ", rank);
    for(int i = 0; i < n*n/(num_proc); i++){
        printf("%lf ", AB[i]);
    }
    printf("\n");
    
    MPI_Type_free(&matrixtype);
    MPI_Finalize();

    return 0;
}