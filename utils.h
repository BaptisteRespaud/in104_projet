#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <math.h>

double calc_mean(double arr[], int n);

double calc_mean_int(int arr[], int n);

double calc_variance(double arr[], int n);

double calc_std(double arr[], int n);  // standard deviation

double *calc_autocorr(double arr[], int n,int n_autoc);  // autocorrelation

void write_to_file_iter(double *arr, int n, char *name);

double *read_from_file(char *name, int n);

double *data_blocking(double *arr, int n, int nblk_max);  // data blocking

void write_blocking_to_file(double *std_blk, int nblk_max, char *name);

#endif