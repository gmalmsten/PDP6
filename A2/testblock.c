#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <cblas.h>
#include <string.h>

#define IND(i,j)  i*n+j
#define N_LIM 10
#define PRODUCE_OUTPUT_FILE

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

int write_output(char *file_name, const double *output, int num_values) {
	FILE *file;
	if (NULL == (file = fopen(file_name, "w"))) {
		perror("Couldn't open output file");
		return -1;
	}
	for (int i = 0; i < num_values; i++) {
		if (0 > fprintf(file, "%.4f ", output[i])) {
			perror("Couldn't write to output file");
		}
	}
	if (0 > fprintf(file, "\n")) {
		perror("Couldn't write to output file");
	}
	if (0 != fclose(file)) {
		perror("Warning: couldn't close output file");
	}
	return 0;
}

void print_matrix(int n,int m, double *A){
    if(n<N_LIM)
    for(int i = 0; i<m; i++){
        for(int j = 0; j<n; j++){
            printf("%lf ", A[i*n+j]);
        }
        printf("\n");
    }
}

void print(int n, double *C){
    for(int i = 0; i < n; i++){
        printf("%lf ", C[i]);
    }
    printf("\n");
}

int main(int argc, char *argv[]){
    if (argc != 1) {
		printf("Usage: %s\n", argv[0]);
		return 1;
	}
	

    // Initialize MPI
	MPI_Init(&argc, &argv);

	int rank, num_proc, left, right;
	MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Create periodic virtual topology
	MPI_Comm GRID_COMM;
	int dims[1];
	int periods[1];
	int reorder = 0;  
	dims[0] = num_proc;
	periods[0] = 1;

	MPI_Cart_create(MPI_COMM_WORLD, 1, dims, periods, reorder, &GRID_COMM);
    MPI_Cart_shift(GRID_COMM, 0, -1, &right, &left);


    // Initialize matrices
    double *A, *B, *C;
    int n = 8, m;
    // Source process reads input
    
    // strart time
    double start_time = MPI_Wtime();

    // Broadcast set up parameters
    m = n/num_proc;
    int avg_square = m*m;
    int rect_size = m*n;
    int matrix_size = n*n;

    printf("m: %d n: %d\n", m, n);
    
    MPI_Datatype rowtype;
    MPI_Type_vector(m, n, n, MPI_DOUBLE, &rowtype);
    MPI_Type_commit(&rowtype);
    
    MPI_Datatype columntype;
    MPI_Type_vector(n, m, n, MPI_DOUBLE, &columntype);
    MPI_Type_commit(&columntype);

    MPI_Datatype squaretype;
    MPI_Type_vector(m, m, n, MPI_DOUBLE, &squaretype);
    MPI_Type_commit(&squaretype);

    MPI_Status status;

    // Allocate memory for local result
    double *C_temp = (double *)malloc(rect_size*sizeof(double));

    for(int i = 0; i < rect_size; i++){
        int offset = (rank==0)*rect_size;
        C_temp[i] = i + offset;
    }

    // Allocate memory for result
    if(rank==0){
        C = (double *)calloc(2*matrix_size, sizeof(double));
    }
    MPI_Request requests[num_proc];
    if(rank == 1){
        for(int i = 0; i<num_proc; i++){
        MPI_Send(&C_temp[i*rect_size/num_proc], rect_size/num_proc, MPI_DOUBLE, 0, i, GRID_COMM);
        }
    }
    printf("Send complete\n");

    // MPI_Type_vector(m, m, n, MPI_DOUBLE, &squaretype);

    if(rank == 0){
        for(int i = 0; i < num_proc; i++){
        MPI_Recv(&C[m*i], 1, squaretype, 1, i, GRID_COMM, &status);
        }
        printf("Received from %d\n", 1);
        
    }

    // Wait for process 0 to receive all local results
    if(rank==0){
        print_matrix(n, n, C);
    }

    MPI_Finalize();

    

    return 0;
}