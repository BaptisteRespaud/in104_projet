#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <math.h>

double calc_mean(double arr[], int n){
    double sum = 0;
    for (int i = 0; i < n; i++){
        sum += arr[i];
    }
    return (double) sum / (double) n;
}

double calc_mean_int(int arr[], int n){
    int sum = 0;
    for (int i = 0; i < n; i++){
        sum += arr[i];
    }
    return (double) sum / (double) n;
}

double calc_variance(double arr[], int n){
    double mean = calc_mean(arr, n);
    double sum = 0.0;
    for (int i = 0; i < n; i++){
        double diff = arr[i] - mean;
        sum += diff * diff;
    }
    return sum / n;
}

double calc_std(double arr[], int n){
    double res = calc_variance(arr, n);
    return sqrt(res);
}

double* calc_autocorr(double arr[], int n, int n_autoc) {
    double* autocorr = malloc(n_autoc * sizeof(double));

    double mean = calc_mean(arr, n);
    double variance = calc_variance(arr, n);  

    if (variance == 0) {
        // specific case autocorr[0] = 1 and forall other index i autocorr[i] = 0
        for (int j = 0; j < n_autoc; j++)
            autocorr[j] = (j == 0) ? 1.0 : 0.0;
        return autocorr;
    }

    for (int j = 0; j < n_autoc; j++) {
        double sum = 0.0;
        for (int i = 0; i < n - j; i++) {
            sum += (arr[i] - mean) * (arr[i + j] - mean);
        }
        autocorr[j] = (double) sum / ((n - j) * variance);
    }

    return autocorr;
}


void write_to_file_iter(double* arr, int n, char* name){
    FILE* f = fopen(name, "w");
    for(int i = 0; i < n; i++){
        fprintf(f, "%lf\n", arr[i]);
    }
    fclose(f);
}

double* read_from_file(char* name, int n){
    double* res = malloc(n * sizeof(double));
    FILE* f = fopen(name, "r");
    double elem;
    for(int i = 0; i < n; i++){
        fscanf(f, "%lf", &elem);
        res[i] = elem;
    }
    fclose(f);
    return res;
}
