#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(int argc, char * argv[]) {
    if(argc!=2){
        return -1;
    }
    int num_integers = atoi(argv[1]);
    char filename[20];
    sprintf(filename, "Input/%d.txt", num_integers);

    srand(time(NULL)); // seed the random number generator

    FILE* fp = fopen(filename, "w");
    fprintf(fp, "%d ", num_integers);
    for (int i = 0; i < num_integers; i++) {
        int rand_int = rand() % (num_integers / 10) + 1;
        fprintf(fp, "%d ", rand_int);
    }
    fclose(fp);

    return 0;
}
