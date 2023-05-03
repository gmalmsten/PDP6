#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

int cmp(const void *a, const void *b)
{
    return (*(int *)a>*(int *)b);
}

int split(int *list, const int pivot){
    int split_index = -1;
    do{
        split_index++;
    } while(list[split_index] < pivot);
    return split_index;
}
void geek_merge(int arr[], int l, int m, int r)
{
    // From https://www.geeksforgeeks.org/merge-sort/
    int i, j, k;
    int n1 = m - l + 1;
    int n2 = r - m;
 
    /* create temp arrays */
    int L[n1], R[n2];
 
    /* Copy data to temp arrays L[] and R[] */
    for (i = 0; i < n1; i++)
        L[i] = arr[l + i];
    for (j = 0; j < n2; j++)
        R[j] = arr[m + 1 + j];
 
    /* Merge the temp arrays back into arr[l..r]*/
    i = 0; // Initial index of first subarray
    j = 0; // Initial index of second subarray
    k = l; // Initial index of merged subarray
    while (i < n1 && j < n2) {
        if (L[i] <= R[j]) {
            arr[k] = L[i];
            i++;
        }
        else {
            arr[k] = R[j];
            j++;
        }
        k++;
    }
 
    /* Copy the remaining elements of L[], if there
    are any */
    while (i < n1) {
        arr[k] = L[i];
        i++;
        k++;
    }
 
    /* Copy the remaining elements of R[], if there
    are any */
    while (j < n2) {
        arr[k] = R[j];
        j++;
        k++;
    }
}


int* merge(int *list1, int *list2, int n1, int n2)
{
    // With inspiration from https://www.geeksforgeeks.org/merge-sort/
    // Initial indeces of list1, list2 and merged list
    int i = 0, j = 0, k = 0;
    int *list = (int *)malloc((n1 + n2)*sizeof(int));
    while(i < n1 && j < n2)
    {
        if(list1[i] <= list2[j]){
            list[k] = list1[i];
            i++;
        }
        else
        {
            list[k] = list2[j];
            j++;
        }
        k++;
    }

    // Copy remaining elements of list1 and list2
        while (i < n1) {
        list[k] = list1[i];
        i++;
        k++;
    }
    while (j < n2) {
        list[k] = list2[j];
        j++;
        k++;
    }

    return list;
}

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
    int n, max_iter, limit;
    max_iter = log2(num_proc);

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
    
    // Divide the data into p ~equal parts 
    MPI_Scatterv(global_list, chunks, displs, MPI_INT, local_list, chunk, MPI_INT, 0, MPI_COMM_WORLD);
    

    // Sort the data locally
    qsort(local_list, chunks[rank], sizeof(int), cmp);

    for(int iter = 0; iter < max_iter; iter++)
    {
        limit = num_proc/2; // rank < limit responsible for list smaller than pivot
        int friend = (rank + limit)%num_proc;
        printf("rank %d friend %d\n", rank, friend);

        
        // Perform global sort
        // Select pivot element within each processor set
        int pivot;
        if(rank ==0){
            // Choose pivot strategy
            pivot = local_list[chunks[0]/2];    // Median of communicator(rank 0)
        }
        MPI_Bcast(&pivot, 1, MPI_INT, 0, MPI_COMM_WORLD);
        int split_index = split(local_list, pivot);


        // Send and receive sub lists
        MPI_Status status;
        int n1 = chunks[rank], n2;             // Sizes of sub lists
        MPI_Sendrecv(&n1, 1, MPI_INT, friend, rank, &n2, 1, MPI_INT, friend, friend, MPI_COMM_WORLD, &status);
        local_list = (int *)realloc(local_list, (n1 + n2)*sizeof(int));
        MPI_Sendrecv(local_list, n1, MPI_INT, friend, rank, &local_list[n1], n2, MPI_INT, friend, friend, MPI_COMM_WORLD, &status);

        // Merge sub lists
        geek_merge(local_list, 0, n1, n1+n2);
    }
    print_list(local_list, chunk, rank);

    MPI_Finalize();

    return 0;
}

// void qsort(void *base, size_t nitems, size_t size, int (*compar)(const void *, const void*))
