#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <cblas.h>

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
	// MPI_Comm GRID_COMM;
	// int dims[2];
	// int periods[2];
	// int reorder = 0;  // May use reorder and MPI_Cart_rank(GRID_COMM, coords[], &rank);
	// dims[0] = sqrt(num_proc);
    // dims[1] = sqrt(num_proc);
	// periods[0] = 0;
    // periods[1] = 0;

	// MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &GRID_COMM);

    // // Check that the processes have the correct coords
    // int coords[2];
    // MPI_Cart_coords(GRID_COMM, rank, 2, coords);
    // printf("Rank %d coords (%d, %d)\n", rank, coords[0], coords[1]);


    // Create periodic virtual topology
	MPI_Comm GRID_COMM;
	int dims[1];
	int periods[1];
	int reorder = 0;  
	dims[0] = num_proc;
	periods[0] = 0;

	MPI_Cart_create(MPI_COMM_WORLD, 1, dims, periods, reorder, &GRID_COMM);


    // Initialize matrices
    double *A, *B, *C, *input;
    int n, m;

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

    // Transpose B matrix to send rows using scatter
    for(int i = 0; i < n; i++)
        for(int j = i + 1; j < n; j++){
            double temp = input[n*n + n*i + j];
            input[n*n + n*i + j] = input[n*n + n*j + i];
            input[n*n + n*j + i] = temp;
        }

    printf("\n");
    printf("B': \n");
    for(int i = 0; i<n; i++){
        for(int j = 0; j<n; j++){
            printf("%lf ", input[n*n + n*i + j]);
        }
        printf("\n");
    }
    
    }

    // Broadcast set up parameters
    MPI_Bcast(&n, 1, MPI_INT, 0, GRID_COMM);
    m = n/num_proc;

    // Define matrix type
    MPI_Datatype matrixtype;
    MPI_Type_vector(m, n, n, MPI_DOUBLE, &matrixtype);

    // MPI_Datatype matrixtype;
    // MPI_Type_vector(sub_n, sub_n, n, MPI_DOUBLE, &matrixtype);
    // MPI_Type_commit(&matrixtype);


    // Scatter matrices to processes
    A = (double *)malloc((m*n)*sizeof(double));
    B = (double *)malloc((n*m)*sizeof(double));
    MPI_Scatter(input, 1, matrixtype, A, n*m, MPI_DOUBLE, 0, GRID_COMM);
    MPI_Scatter(&input[n*n], 1, matrixtype, B, n*m, MPI_DOUBLE, 0, GRID_COMM);


    // // Print submatrices as rows


    // Allocate memory for local result
    C = (double *)malloc(m*m*sizeof(double));
    if(rank==0){
        printf("m: %d n: %d\n", m, n);
    for(int i = 0; i < m; i++){
    for(int j = 0; j < n; j++){
        printf("%lf ", A[i*n + j]);
    }
    printf("\n");
    }
    printf("\n");

    for(int i = 0; i < n; i++){
    for(int j = 0; j<m; j++){
        printf("%lf ", B[i + j*n]);
    }
    printf("\n");
    }
    printf("\n");
        double *CTEST = (double *)malloc(m*m*sizeof(double));
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, m, n, 1.0, A, n, B, m, 0.0, CTEST, m);
        for(int i = 0; i < m; i++){
            for(int j = 0; j < m; j++){
                printf("%lf ",CTEST[i*m + j]);
            }
            printf("\n");
        }
        
        free(CTEST);
    }
    
    MPI_Type_free(&matrixtype);
    MPI_Finalize();

    return 0;
}