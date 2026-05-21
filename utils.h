#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <math.h>

// compute the mean of the double array arr
double calc_mean(double arr[], int n);

// compute the mean of the double array arr
double calc_mean_int(int arr[], int n);

// compute the variance of the double array arr
double calc_variance(double arr[], int n);

// compute the standard deviation of the double array arr
double calc_std(double arr[], int n); 

// compute the autocorrelation array of size n_autoc of a double array arr of size n
double *calc_autocorr(double arr[], int n,int n_autoc);  // autocorrelation

/* write down in file named name the array arr following this syntaxe 
arr[0]
arr[1]
...
arr[n-1]
*/
void write_to_file_iter(double *arr, int n, char *name);

/* produce an array res of size n from the file named name following this syntaxe 
arr[0]
arr[1]
...
arr[n-1] 
*/
double *read_from_file(char *name, int n);

#endif