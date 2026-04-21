#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <math.h>

// compute the mean of the double array arr
double calc_mean(double arr[], int n){
    double sum = 0;
    for (int i = 0; i < n; i++){
        sum += arr[i];
    }
    return sum/n;
}

// compute the mean of the int array arr
double calc_mean_int(int arr[], int n){
    int sum = 0;
    for (int i = 0; i < n; i++){
        sum += arr[i];
    }
    return sum / n;
}

// compute the variance of the double array arr
double calc_variance(double arr[], int n){
    double mean = calc_mean(arr, n);
    int sum = 0;
    for (int i = 0; i < n; i++){
        sum += (arr[i] - mean) * (arr[i] - mean);
    }
    return sum/n;
}

// compute the standard deviation of the double array arr
double calc_std(double arr[], int n){
    double res = calc_variance(arr, n);
    return sqrt(res);
}

/* write down in file named name the array arr following this syntaxe 
arr[0]
arr[1]
...
arr[n-1]
*/
void write_to_file_iter(double* arr, int n, char* name){
    FILE* f = fopen(name, "w");
    for(int i = 0; i < n; i++){
        fprintf(f, "%lf\n", arr[i]);
    }
    fclose(f);
}

/* produce an array res of size n from the file named name following this syntaxe 
arr[0]
arr[1]
...
arr[n-1] 
*/
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

double *calc_autocorr(double arr[], int n,int n_autoc){
    double* autocorr_arr = malloc(n_autoc * sizeof(double));
    double mean = calc_mean(arr, n);
    double variance = calc_variance(arr, n);
    for (int j = 0; j < n_autoc; j++){
        double* prod = malloc((n-j) * sizeof(double));
        double retard = malloc((n-j) * sizeof(double));
        for (int i = 0; i <n-n_autoc; i++){
            prod[i] = arr[i]*arr[i+j];
            retard[i]= arr[i+j];
        }
        autocorr_arr[j] = (calc_mean(prod,n-j) - calc_mean(retard,n-j)*mean/variance);
        free(prod);
        free(retard);
    }
    return autocorr_arr;