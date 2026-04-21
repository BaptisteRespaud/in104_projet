#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include "../utils.h"

#define PI 3.1415926535

int main(){

    // test of read_from_file and write_to_file_iter
    char* name = "test_read_write_from_file";
    int n = 1000;
    
    double* res = malloc(n*sizeof(double));

    for (int i = 0; i < n; i++){
        res[i] = (double) rand() / (double) RAND_MAX;
    }

    write_to_file_iter(res, n, name);

    double* read = read_from_file(name, n);

    if (n != 0){
        for(int i = 0; i < n; i++){
            assert(res[i] - read[i] <= 1);
        }
    }

    // test of the autocorrellation function 
    int N = 100;
    double dt = 4 * PI / N;

    double* x = malloc(N * sizeof(double));

    for (int i = 0; i < N; i++) {
        double noise = 0.1 * (2.0 * rand() / RAND_MAX - 1.0); // noise ±0.1
        x[i] = sin(dt * i) + noise;
    }

    int n_autoc = N / 2;

    double* ac = calc_autocorr(x, N, n_autoc);

    for (int j = 0; j < n_autoc; j++)
        printf("ac[%d] = %lf\n", j, ac[j]);

    free(ac);
    free(x);

}