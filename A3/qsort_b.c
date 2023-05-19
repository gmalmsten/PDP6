#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

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

int cmp(const void *a, const void *b)
{
    return (*(int *)a > *(int *)b);
}

// Funciton that finds the split index of a sorted list
int split(int *list, const int pivot, const int n)
{
    int split_index = -1;
    do
    {
        split_index++;
    } while ((split_index < n) && (list[split_index] < pivot));
    return split_index;
}

void merge(int *list1, int *list2, int n1, int n2, int *list)
{
    // With inspiration from https://www.geeksforgeeks.org/merge-sort/
    // Initial indeces of list1, list2 and merged list
    int i = 0, j = 0, k = 0;

    while (i < n1 && j < n2 )
    {
        if (list1[i] <= list2[j])
        {
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

    //Copy remaining elements of list1 and list2
    while (i < n1)
    {
        list[k] = list1[i];
        i++;
        k++;
    }
    while (j < n2)
    {
        list[k] = list2[j];
        j++;
        k++;
    }
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
            printf("i: %d\n", i);
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

int write_output(char *file_name, const int *output, int num_values) {
	FILE *file;
	if (NULL == (file = fopen(file_name, "w"))) {
		perror("Couldn't open output file");
		return -1;
	}
	for (int i = 0; i < num_values; i++) {
		if (0 > fprintf(file, "%d ", output[i])) {
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

int main(int argc, char *argv[])
{
    if (argc != 4)
    {
        printf("Usage %s n output_file pivot_strategy\n", argv[0]);
        return -1;
    }

    const int n = atoi(argv[1]);
    char *output_file = argv[2];
    int pivot_strategy = atoi(argv[3]);

    // Parallel environment
    MPI_Init(&argc, &argv);
    int global_rank, global_num_proc;

    MPI_Comm_rank(MPI_COMM_WORLD, &global_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &global_num_proc);

    int *global_list, *local_list;
    int max_iter, limit;
    max_iter = log2(global_num_proc);

    if (global_rank == 0)
    {
        if (n < 0)
        {
            printf("Error n < 0\n");
            return -1;
        }
        printf("%d\t%d\t", pivot_strategy, n);
        double start = MPI_Wtime();
        global_list = (int *)malloc(n*sizeof(int));
        for(int i = 0; i < n; i++){
            global_list[i] = n-i;
        }
        
    }


    // Sizes of chunks in processors
    int chunk = n / global_num_proc;     
    int remainder = n % global_num_proc; 

    // Allocate memory for storing chunks and displacements, for Scatterv/Gatherv
    int *chunks = malloc(sizeof(int) * global_num_proc);
    int *displs = malloc(sizeof(int) * global_num_proc);

    // Calculate chunks and displacements
    for (int i = 0; i < global_num_proc; i++)
    {
        chunks[i] = chunk;
        displs[i] = i * chunk;
    }

    // Add remainder to last chunk
    if (global_rank == global_num_proc - 1)
    {
        chunk += remainder;
        local_list = malloc(sizeof(int) * chunk); // Allocate memory for last chunk
    }
    else
        local_list = malloc(sizeof(int) * chunk); // Allocate memory for all other chunks

    chunks[global_num_proc - 1] += remainder; // Add remainder to last chunk

    

    // Create new communicator to split in iterations
    MPI_Comm MPI_LOCAL_COMM;
    MPI_Comm_split(MPI_COMM_WORLD, 0, global_rank, &MPI_LOCAL_COMM);
    int rank, num_proc;
    
    chunk = chunks[global_rank];
    int *send_list = (int *)malloc(chunk* sizeof(int));
    int *remaining_list = (int *)malloc(chunk * sizeof(int));
    

    // --------- START TIMER --------- //
    double start_time = MPI_Wtime();
    // --------- START TIMER --------- //

    // Divide the data into p ~equal parts
    MPI_Scatterv(global_list, chunks, displs, MPI_INT, local_list, chunk, MPI_INT, 0, MPI_COMM_WORLD);

    // Sort the data locally
    qsort(local_list, chunks[global_rank], sizeof(int), cmp);

    // Perform global sort
    for (int iter = 0; iter < max_iter; iter++)
    {
        MPI_Comm_rank(MPI_LOCAL_COMM, &rank);
        MPI_Comm_size(MPI_LOCAL_COMM, &num_proc);
        
        limit = num_proc / 2; // rank < limit responsible for list smaller than pivot
        int friend = (rank + limit) % num_proc;
        
        // Select pivot element within each processor set
        int pivot, median;
        switch (pivot_strategy) {

            case 1:
                // Select the median in one processor in each group of processors.
                if (rank == 0)
                {
                    pivot = local_list[chunk / 2];
                }
                
                break;

            case 2:
                //Select the median of all medians in each processor group.
                median = local_list[chunk/2];
                int *medians = (int *)malloc(num_proc*sizeof(int));
                MPI_Gather(&median, 1, MPI_INT, medians, 1, MPI_INT, 0, MPI_LOCAL_COMM);
                if(rank == 0)
                {
                    qsort(medians, num_proc, sizeof(int), cmp);
                    pivot = medians[num_proc/2];  
                }
                break;

            case 3:
                //Select the mean value of all medians in each processor group.
                median = local_list[chunk/2];
                int sum;
                MPI_Reduce(&median, &sum, 1, MPI_INT, MPI_SUM, 0, MPI_LOCAL_COMM);
                pivot = sum / num_proc;
                break;
                
            default:
                if(rank == 0)
                {
                    printf("Invalid pivot strategy\n");
                    printf("Pivot strategy options:\n 1 - Median of middle communicator \n 2 - Median of medians\n 3 - Mean of medians\n");
                    return -1;
                }

        }
        MPI_Bcast(&pivot, 1, MPI_INT, 0, MPI_LOCAL_COMM);
        int split_index = split(local_list, pivot, chunk);

        // Send and receive sub lists
        MPI_Status status;
        int send_n, receive_n; 
        if (rank < limit) // Send upper part of sub array if rank is in lower part
        {
            // Sending upper part of sub array
            send_n = chunk - split_index;

            send_list = (int *)realloc(send_list, send_n*sizeof(int));
            remaining_list = (int *)realloc(remaining_list, split_index*sizeof(int));
            
            memcpy(send_list, &local_list[split_index], send_n * sizeof(int));
            memcpy(remaining_list, local_list, split_index * sizeof(int));
        }
        else
        {
            // Sending lower part of sub array
            send_n = split_index;
            send_list = (int *)realloc(send_list, send_n*sizeof(int));
            remaining_list = (int *)realloc(remaining_list, (chunk - split_index)*sizeof(int));

            memcpy(send_list, local_list, send_n * sizeof(int));
            memcpy(remaining_list, &local_list[split_index], (chunk - split_index) * sizeof(int));
        }
        int remaining_chunk = chunk - send_n;

        // Sendrecv sizes of subarrays to recv
        MPI_Sendrecv(&send_n, 1, MPI_INT, friend, rank, &receive_n, 1, MPI_INT, friend, friend, MPI_LOCAL_COMM, &status);
        int *received_list = (int *)malloc(receive_n * sizeof(int));

        // Sendrecv subarrays
        MPI_Request request;
        MPI_Isend(send_list, send_n, MPI_INT, friend, rank, MPI_LOCAL_COMM, &request);
        MPI_Recv(received_list, receive_n, MPI_INT, friend, friend, MPI_LOCAL_COMM, &status);
        
        // Update size of local array
        chunk = chunk - send_n + receive_n;
        printf("Rank %d chunk %d\n", rank, chunk);
        if(chunk > 0){
            local_list = (int *)realloc(local_list, (chunk) * sizeof(int));
            merge(remaining_list, received_list, remaining_chunk, receive_n, local_list);
        }
        
        // Split the communicator into two halves, those with sub array
        // below pivot and those with values above pivot
        MPI_Comm_split(MPI_LOCAL_COMM, rank < limit, rank%limit, &MPI_LOCAL_COMM);
        MPI_Comm_rank(MPI_LOCAL_COMM, &rank);
        
    }


    // Clean up
    free(send_list); 
    free(remaining_list);

    // Gather the locally sorted lists
    MPI_Gather(&chunk, 1, MPI_INT, chunks, 1, MPI_INT, 0, MPI_COMM_WORLD);
    for(int p = 1; p < global_num_proc; p++){
        displs[p] = displs[p-1] + chunks[p-1];
    }
    MPI_Gatherv(local_list, chunk, MPI_INT, global_list, chunks, displs, MPI_INT, 0, MPI_COMM_WORLD);
    
    // --------- END TIMER --------- //
    double end_time = MPI_Wtime();
    double time = end_time - start_time;
    double max_time;
    MPI_Reduce(&time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (global_rank == 0)
    {
        printf("%f\n", max_time);
    }
    // --------- END TIMER --------- //


    // Check results and print to output file
    if(global_rank == 0){
        for(int i = 0; i < n - 1; i++)
        {
            if(global_list[i] > global_list[i+1]){
                printf("List not sortered\n");
                return -1;
            }
        }
        #ifdef PRODUCE_OUTPUT_FILE
            write_output(output_file, global_list, n);
        #endif
    }


    // free memory
    free(chunks);
    free(displs);
    free(global_list);
    free(local_list);

    // free communicator
    MPI_Comm_free(&MPI_LOCAL_COMM); 

    MPI_Finalize();

    return 0;
}
