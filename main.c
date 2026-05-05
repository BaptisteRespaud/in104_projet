#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <math.h>
#include "utils.h"
#include "lattice.h"
#include "model.h"
#include "metropolis.h"

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
        if (r < 0.5){
            state[i] = -1;
        } else {
            state[i] = 1; 
        }
    }

    // metropolis algorithm 

    int T_iter = DeltaT/dT;

    double* sqr_avg_mag_vsT = malloc(T_iter * sizeof(double));
    double* sqr_avg_mag_arr = malloc(N_tot * sizeof(double));


    for (int j = 0 ; j < T_iter; j ++){
        // equilibration steps
        for (int k = 0; k < N_eq; k++){
            move(model, state);    
        } 

        // computing of the observable
        model->T = T + j * dT;
        for (int i = 0 ; i < N_tot; i++){

            // compute the average squared magnetization 
            double m = avg_magnetization(model, state);
            sqr_avg_mag_arr[i] = m * m;
            move(model, state);
        }

        sqr_avg_mag_vsT[j] = calc_mean(sqr_avg_mag_arr, N_tot);
        // printf("average squared magnetization after %d iterations and with T = %f: %lf\n", N_tot, T + j*dT, sqr_avg_mag_vsT[j]);
    }

    char filename[64];
    sprintf(filename, "iterations/sqr_avg_mag_vsT_L%d.txt", Lx);
    write_to_file_iter(sqr_avg_mag_vsT, T_iter, filename);

    free(model->lat);
    free(model);
    free(state);
    free(sqr_avg_mag_arr);
    free(sqr_avg_mag_vsT);
    return 0;
}


