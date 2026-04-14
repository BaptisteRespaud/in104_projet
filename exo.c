#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include "utils.h"

int main(){
    // estimation of Pi

    srand(time(NULL));

    int N = 100000;
    int in_the_circle = 0;
    double* pi_arr = malloc(N*sizeof(double));

    for (int i = 0; i < N; i++){
        double x = (double) rand() / (double) RAND_MAX;
        double y = (double) rand() / (double) RAND_MAX;
        if (x * x + y * y <= 1){
            in_the_circle += 1;
        }
        pi_arr[i] = 4 * ((double) in_the_circle / (double) N);
    }

    printf("pi = %f\n", pi_arr[N-1]);
    write_to_file_iter(pi_arr, N, "pi_iterations");
    
    // estimation of the Gaussian integral

    N = 10000;
    int count = 0;

    for (int i = 0; i < N; i++){
        double x = 8.0 * rand() / RAND_MAX - 4.0;

        double y = (double) rand() / (double) RAND_MAX;
        if (y <= exp(- x*x / 2)){
            count += 1;
        }
    }

    double estimation_of_gaussian = 8 * (double) count / (double) N;
    printf("gaussian = %f\n", estimation_of_gaussian);
    return 0;
}