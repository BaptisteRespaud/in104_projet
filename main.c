#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <math.h>
#include "utils.h"
#include "lattice.h"
#include "model.h"
#include "metropolis.h"
#include "swendsen_wang.h"

int main() {
    srand(time(NULL));

    // read params from file "config.txt"

    double T, B;
    int Lx, Ly, N_tot, N_eq;
    double DeltaT, dT;

    FILE *fp = fopen("config.txt", "r");
    fscanf(fp, "T %lf \n", &T);
    printf("Initial Temperature T = %lf \n", T);
    
    fscanf(fp, "B %lf \n", &B);
    printf("External field B = %lf \n", B);
    
    fscanf(fp, "Lx %d \n", &Lx);
    fscanf(fp, "Ly %d \n", &Ly);
    printf("Square lattice Lx = %d Ly = %d \n", Lx, Ly);
    
    fscanf(fp, "N_tot %d \n", &N_tot);
    fscanf(fp, "N_eq %d \n", &N_eq);
    printf("Metropolis N_tot = %d \n", N_tot);
    printf("Metropolis N_eq = %d \n", N_eq);

    fscanf(fp, "DeltaT %lf \n", &DeltaT);
    fscanf(fp, "dT %lf \n", &dT);
    printf("Range of temperatures DeltaT = %lf\n", DeltaT);
    printf("Integration step dT = %lf\n", dT);
    fclose(fp);

    // initialization 

    SquareLattice* lat = malloc(sizeof(SquareLattice));
    lat->L = Lx*Ly; lat->Lx = Lx; lat->Ly = Ly; lat->pbc = true;
    ising* model = define_model(1, B, T, lat);

    int* state = malloc(lat->L * sizeof(int));
    for(int i = 0; i < lat->L; i++){
        double r = (double) rand() / (double) RAND_MAX;
        state[i] = (r < 0.5 ? -1 : 1);
    }

    // metropolis algorithm 

    int T_iter = DeltaT/dT;
    int N = Lx * Ly;   // number of spins

    double* sqr_avg_mag_vsT = malloc(T_iter * sizeof(double));
    double* sqr_avg_mag_arr = malloc(N_tot * sizeof(double));
    double* m4_arr = malloc(N_tot * sizeof(double));

    double* E_arr = malloc(N_tot * sizeof(double));
    double* E_sqr_arr = malloc(N_tot * sizeof(double));
    double* avg_E_vsT = malloc(T_iter * sizeof(double));
    double* cv_vsT = malloc(T_iter* sizeof(double));
    double* binder_vsT = malloc(T_iter * sizeof(double));

    for (int j = 0 ; j < T_iter; j ++){
        model->T = T + j * dT;

        // equilibration steps
        /*for (int k = 0; k < N_eq; k++){
            //---metropolis---
            move(model, state);    
        } 
        */

        // computing of the observables
        for (int i = 0 ; i < N_tot; i++){
            double m = avg_magnetization(model, state);
            double m2 = m * m;
            double m4 = m2 * m2;

            double E = total_energy(model, state);

            sqr_avg_mag_arr[i] = m2;
            m4_arr[i] = m4;

            E_arr[i] = E;
            E_sqr_arr[i] = E * E;

            //---metropolis---
            //move(model, state);

            //---swendsen_wang---
            swendsen_wang(model, state);
        }

        double avg_m2 = calc_mean(sqr_avg_mag_arr, N_tot);
        double avg_m4 = calc_mean(m4_arr, N_tot);

        sqr_avg_mag_vsT[j] = avg_m2;

        double avg_E = calc_mean(E_arr, N_tot);
        double avg_E_sqr = calc_mean(E_sqr_arr, N_tot);

        avg_E_vsT[j] = avg_E;

        // Correct heat capacity per spin
        cv_vsT[j] = (avg_E_sqr - avg_E * avg_E) / (N * model->T * model->T);

        // Correct Binder cumulant of magnetization
        binder_vsT[j] = 1 - avg_m4 / (3 * avg_m2 * avg_m2);
    }
    
    char filename1[64];
    sprintf(filename1, "iterations/sqr_avg_mag_vsT_L%d.txt", Lx);
    write_to_file_iter(sqr_avg_mag_vsT, T_iter, filename1);

    char filename2[64];
    sprintf(filename2, "iterations/avg_E_vsT_L%d.txt", Lx);
    write_to_file_iter(avg_E_vsT, T_iter, filename2);

    char filename3[64];
    sprintf(filename3, "iterations/Cv_vsT_L%d.txt", Lx);
    write_to_file_iter(cv_vsT, T_iter, filename3);

    char filename4[64];
    sprintf(filename4, "iterations/binder_vsT_L%d.txt", Lx);
    write_to_file_iter(binder_vsT, T_iter, filename4);

    // compute the autocorrelation function for T = 2.3J
    model->T = 2.3;

    for (int k = 0; k < N_eq; k++){
        move(model, state);    
    } 

    for (int i = 0 ; i < N_tot; i++){
        double m = avg_magnetization(model, state);
        sqr_avg_mag_arr[i] = m * m;
        move(model, state);
    }

    double* autocorr_23_mag = calc_autocorr(sqr_avg_mag_arr, N_tot, 1000);
    char filename5[64];
    sprintf(filename5, "iterations/autocorr_mag_L%d.txt", Lx);
    write_to_file_iter(autocorr_23_mag, 1000, filename5);

    free(model->lat);
    free(model);
    free(state);
    free(sqr_avg_mag_arr);
    free(sqr_avg_mag_vsT);
    free(m4_arr);
    free(E_arr);
    free(E_sqr_arr);
    free(cv_vsT);
    free(avg_E_vsT);
    free(autocorr_23_mag);
    free(binder_vsT);

    return 0;
}
