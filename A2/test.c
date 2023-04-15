#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

int main(int argc, char** argv) {
    int rank, size;
    int n = 4; // size of the matrix
    int p = 4; // number of elements
    int sub_n; // size of the sub-matrix for each process
    int *matrix = NULL; // the matrix
    int *sub_matrix = NULL; // sub-matrix for each process
    int *counts = NULL; // number of elements to send to each process
    int *displs = NULL; // displacement for each process in the matrix

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // make sure that p divides n^2
    if (n * n % p != 0) {
        printf("Error: p must divide n^2.\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // determine the size of the sub-matrix for each process
    sub_n = n / sqrt(p);

    // allocate memory for the matrix and sub-matrix
    if (rank == 0) {
        matrix = (int*) malloc(n * n * sizeof(int));
        for (int i = 0; i < n * n; i++) {
            matrix[i] = i;
        }
    }
    sub_matrix = (int*) malloc(sub_n * sub_n * sizeof(int));

    // calculate the counts and displacements for each process
    counts = (int*) malloc(size * sizeof(int));
    displs = (int*) malloc(size * sizeof(int));
    for (int i = 0; i < size; i++) {
        counts[i] = sub_n * sub_n;
        displs[i] = (i / sqrt(p)) * sub_n * n + (i % (int)sqrt(p)) * sub_n;
    }

    // scatter the sub-matrix to all processes
    MPI_Scatterv(matrix, counts, displs, MPI_INT, sub_matrix, sub_n * sub_n, MPI_INT, 0, MPI_COMM_WORLD);

    // print the sub-matrix for each process
    printf("Sub-matrix in process %d:\n", rank);
    for (int i = 0; i < sub_n; i++) {
        for (int j = 0; j < sub_n; j++) {
            printf("%d ", sub_matrix[i * sub_n + j]);
        }
        printf("\n");
    }

    // free memory and finalize MPI
    if (rank == 0) {
        free(matrix);
    }
    free(sub_matrix);
    free(counts);
    free(displs);
    MPI_Finalize();
    return 0;
}
