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

double calc_variance(double arr[], int n){
    double x_barre= calc_mean(arr,n);
    int sum=0;
    for (int i=0; i<n; i++){
        sum += (arr[i]-x_barre)**2;
    }
    return sum;
}