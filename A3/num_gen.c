#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main() {
    int num_integers = 1000000;
    char filename[20];
    sprintf(filename, "%d.txt", num_integers);

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
