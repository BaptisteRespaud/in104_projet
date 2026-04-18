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
    return sum/n;
}

double calc_mean_int(int arr[], int n){
    int sum = 0;
    for (int i = 0; i < n; i++){
        sum += arr[i];
    }
    return sum / n;
}

double calc_variance(double arr[], int n){
    double mean = calc_mean(arr, n);
    int sum = 0;
    for (int i = 0; i < n; i++){
        sum += (arr[i] - mean) * (arr[i] - mean);
    }
    return sum/n;
}

double calc_std(double arr[], int n){
    double res = calc_variance(arr, n);
    return sqrt(res);
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
