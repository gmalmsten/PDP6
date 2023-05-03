#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int print_list(int *list, int n, int rank)
{
    printf("Rank %d: ", rank);
    for (int i = 0; i < n; i++)
    {
        printf("%d ", list[i]);
    }
    printf("\n");
    return 0;
}

int read_input(const char *file_name, int **input)
{
    FILE *file;
    if (NULL == (file = fopen(file_name, "r")))
    {
        perror("Couldn't open input file");
        return -1;
    }
    int n;
    if (EOF == fscanf(file, "%d", &n))
    {
        perror("Couldn't read element count from input file");
        return -1;
    }
    if (NULL == (*input = malloc(n * sizeof(int))))
    {
        perror("Couldn't allocate memory for matrix A");
        return -1;
    }
    for (int i = 0; i < n; i++)
    {
        if (EOF == fscanf(file, "%d", &((*input)[i])))
        {
            perror("Couldn't read elements from input file");
            return -1;
        }
    }
    if (0 != fclose(file))
    {
        perror("Warning: couldn't close input file");
    }
    return n;
}

int main(int argc, char *argv[])
{
    if (argc != 4)
    {
        printf("Usage %s input_file output_file pivot_strategy\n", argv[0]);
        //----------------------------------------------------------------------------------------------------------------------------//
    }

    char *input_file = argv[1];

    // Parallel environment
    MPI_Init(&argc, &argv);
    int rank, num_proc;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);

    int *global_list, *local_list;
    int n;

    if (rank == 0)
    {
        n = read_input(input_file, &global_list);
        if (n == -1)
        {
            printf("Error reading input file\n");
            return -1;
        }
        print_list(global_list, n, rank);
    }

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int chunk = n / num_proc;
    int remainder = n % num_proc;
    int *chunks = malloc(sizeof(int) * num_proc);
    int *displs = malloc(sizeof(int) * num_proc);
    for (int i = 0; i < num_proc; i++)
    {
        chunks[i] = chunk;
        displs[i] = i * chunk;
    }


    if (rank == num_proc - 1)
    {
        chunk += remainder;
        local_list = malloc(sizeof(int) * chunk);
    }
    else
        local_list = malloc(sizeof(int) * chunk);

    chunks[num_proc - 1] += remainder;
    
    MPI_Scatterv(global_list, chunks, displs, MPI_INT, local_list, chunk, MPI_INT, 0, MPI_COMM_WORLD);

    print_list(local_list, chunk, rank);

    MPI_Finalize();

    return 0;
}

// void qsort(void *base, size_t nitems, size_t size, int (*compar)(const void *, const void*))
