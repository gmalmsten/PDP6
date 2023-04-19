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
		if (0 > fprintf(file, "%.6f ", output[i])) {
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
    for(int i = 0; i<n; i++){
        for(int j = 0; j<m; j++){
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
    if (argc != 3) {
		printf("Usage: %s input_file output_file\n", argv[0]);
		return 1;
	}
	char *input_name = argv[1];
	char *output_name = argv[2];

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
    double *A, *B, *C, *input;
    int n, m;
    // Source process reads input
    if(rank == 0){
        n = read_input(input_name, &input);
    }
    

    // Broadcast set up parameters
    MPI_Bcast(&n, 1, MPI_INT, 0, GRID_COMM);
    m = n/num_proc;
    int avg_square = m*m;
    int rect_size = m*n;
    int matrix_size = n*n;
    
    

    MPI_Datatype rowtype;
    MPI_Type_vector(m, n, n, MPI_DOUBLE, &rowtype);
    MPI_Type_commit(&rowtype);
    
    MPI_Datatype columntype;
    MPI_Type_vector(n, m, n, MPI_DOUBLE, &columntype);
    MPI_Type_commit(&columntype);

    MPI_Datatype squaretype;
    MPI_Type_vector(m, m, n, MPI_DOUBLE, &squaretype);
    MPI_Type_commit(&squaretype);


    // Scatter matrices to processes
    A = (double *)malloc((rect_size)*sizeof(double));
    B = (double *)malloc((rect_size)*sizeof(double));
    MPI_Scatter(input, 1, rowtype, A, rect_size, MPI_DOUBLE, 0, GRID_COMM);
    
    // Scatter B matrices, column blocks
    MPI_Request requests[num_proc*num_proc];
    if(rank==0){
        for(int p = 0; p < num_proc; p++){
            MPI_Isend(&input[n*n+p*m], 1, columntype, p, p, GRID_COMM, &requests[p]);
        }
    }
    MPI_Status status;
    MPI_Recv(B, m*n, MPI_DOUBLE, 0, rank, GRID_COMM, &status);
    
    

    // Allocate memory for local result
    double *C_temp = (double *)malloc(rect_size*sizeof(double));


    // Allocate memory for result
    if(rank==0){
        C = (double *)malloc(matrix_size*sizeof(double));
    }

    // Start timer
    double start_time = MPI_Wtime();

    // Do calculations on row and column in memory, shift columns left and repeat
    // Store local results in C_temp, each column row major
    for(int i = 0; i<num_proc; i++)
    {   
        int col = (i+rank)%num_proc;
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, m, n, 1.0, A, n, B, m, 0.0, &C_temp[col*avg_square], m);
        MPI_Sendrecv_replace(B, rect_size, MPI_DOUBLE, left, 0, right, 0, GRID_COMM, &status);
    }

   
    
    // Send local result to root process
    for(int i = 0; i<num_proc; i++){
        int col = (rank + i)%num_proc;
        MPI_Isend(&C_temp[col*avg_square], avg_square, MPI_DOUBLE, 0, rank*num_proc+col, GRID_COMM, &requests[rank*num_proc+col]);
    }

     
    // Gather local results
    if(rank == 0){
        for(int p = 0; p < num_proc; p++){
            for(int col = 0; col<num_proc; col++){
                MPI_Irecv(&C[p*rect_size + m*col], 1, squaretype, p, p*num_proc+col, GRID_COMM, &requests[p*num_proc+col]);
            }
        }
    }

    
    // Wait for process 0 to receive all local results
    if(rank==0){
        MPI_Waitall(num_proc*num_proc, requests, MPI_STATUSES_IGNORE);
        // print_matrix(n,n,C);
    }
     
    
    // Stop timer
    double end_time = MPI_Wtime();
    double time = end_time - start_time;
    double max_time;
    MPI_Reduce(&time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, GRID_COMM);

     
    // Print results
    if(rank == 0){
        printf("%lf\n", max_time);
		#ifdef PRODUCE_OUTPUT_FILE
			if (0 != write_output(output_name, C, n*n)) {
				return 2;
			}
		#endif
        
        // Clean up
        free(C);
        free(input); 
    }
    
    /*
    // Clean up
    free(A);
    free(B);
    free(C_temp);   

    MPI_Type_free(&rowtype);
    MPI_Type_free(&columntype);
    MPI_Type_free(&squaretype);
    MPI_Comm_free(&GRID_COMM);
    */
    MPI_Finalize();

    

    return 0;
}