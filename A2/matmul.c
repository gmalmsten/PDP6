#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <cblas.h>

#define IND(i,j)  i*n+j

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

void print_matrix(int n,int m, double *A){
    for(int i = 0; i<n; i++){
        for(int j = 0; j<m; j++){
            printf("%lf ", A[i*n+j]);
        }
        printf("\n");
    }
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
        
        // // Print read input
        // printf("A: \n");
        // print_matrix(n,n, input);

        // printf("\n");
        // printf("B: \n");
        // for(int i = 0; i<n; i++){
        //     for(int j = 0; j<n; j++){
        //         printf("%lf ", input[n*n + n*i + j]);
        //     }
        //     printf("\n");
        // }
    }

    // strart time
    double start_time = MPI_Wtime();

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
    // MPI_Scatter(&input[n*n], 1, matrixtype, B, n*m, MPI_DOUBLE, 0, GRID_COMM);
    if(rank==0)
    for(int p = 0; p < num_proc; p++){
        MPI_Send(&input[n*n+p*m], 1, columntype, p, p, GRID_COMM);
    }
    MPI_Status status;
    MPI_Recv(B, m*n, MPI_DOUBLE, 0, rank, GRID_COMM, &status);
    free(input);


    // Allocate memory for local result
    C = (double *)malloc(rect_size*sizeof(double));
    double *C_temp = (double *)malloc(avg_square*sizeof(double));

    MPI_Request requests[num_proc*num_proc];

    for(int i = 0; i<num_proc; i++)
    {
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, m, n, 1.0, A, n, B, m, 0.0, C_temp, m);
    // cblas_dgemm (const CBLayout,transa, transb, m, n, k, alpha, *a,lda,*b, ldb, beta,*c, ldc);
    // printf("Rank %d left: %d right %d\n", rank, left, right);
    if(rank == 2){
        print_matrix(m,m, C_temp);
        printf("\n");
        // printf("Process %d sending with tag %d \n", rank, (rank+i)%num_proc);
    }
    // Send local result (avg_square) to process 0, non blocking
    MPI_Isend(C_temp, avg_square, MPI_DOUBLE, 0, rank*num_proc+(i+rank)%num_proc, GRID_COMM, &requests[rank*num_proc+i]);
    // Send and receive columns
    MPI_Sendrecv_replace(B, rect_size, MPI_DOUBLE, left, 0, right, 0, GRID_COMM, &status);
    }


    // Receive all local results in process 0
    if(rank == 0)
    {   
        C = (double *)malloc(matrix_size*sizeof(double));
        for(int p = 0; p < num_proc; p++){
            for(int col = 0; col < num_proc; col++)
            {
                // MPI_Irecv(&C[rect_size*p+col*m], 1, squaretype, p, col, GRID_COMM, &requests[p*num_proc+col]);
                MPI_Status recvstatus;
                MPI_Recv(&C[rect_size*p+col*m], 1, squaretype, p,  p*num_proc+col, GRID_COMM, &recvstatus);
                // printf("RECIEVED %d\n", p);
            }
        }
        // MPI_Waitall(num_proc*num_proc, requests, MPI_STATUSES_IGNORE);
        printf("C: \n");
        print_matrix(n,n,C);
    }



    // end time
    double end_time = MPI_Wtime();
    double time = end_time - start_time;
    double max_time;
    MPI_Reduce(&time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, GRID_COMM);
    if(rank == 0){
        printf("%lf\n", max_time);
    }
    
    // if(rank == 0){
    //     printf("C: ");
    //     print_matrix(n, C);
    // }
    MPI_Type_free(&rowtype);
    MPI_Type_free(&columntype);
    MPI_Type_free(&squaretype);
    MPI_Comm_free(&GRID_COMM);
    MPI_Finalize();

    return 0;
}