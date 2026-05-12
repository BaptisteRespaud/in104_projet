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
    int *count = malloc(sizeof(int));
    *count = 0;

    double* sqr_avg_mag_vsT = malloc(T_iter * sizeof(double));
    double* sqr_avg_mag_arr = malloc(N_tot * sizeof(double));
    double* avg_E_arr = malloc(N_tot * sizeof(double));
    double* avg_E_vsT = malloc(T_iter * sizeof(double));



    for (int j = 0 ; j < T_iter; j ++){
        model->T = T + j * dT;

        // equilibration steps
        for (int k = 0; k < N_eq; k++){
            move(model, state, count);    
        } 

        // computing of the observable
        for (int i = 0 ; i < N_tot; i++){

            // compute the average squared magnetization 
            double m = avg_magnetization(model, state);
            double E = total_energy(model, state);
            sqr_avg_mag_arr[i] = m * m;
            avg_E_arr[i] = E;
            move(model, state, count);
        }
        write_to_file_iter(sqr_avg_mag_arr, N_tot, "test");
        sqr_avg_mag_vsT[j] = calc_mean(sqr_avg_mag_arr, N_tot);
        avg_E_vsT[j] = calc_mean(avg_E_arr, N_tot);
        // printf("average squared magnetization after %d iterations and with T = %f: %lf\n", N_tot, T + j*dT, sqr_avg_mag_vsT[j]);
    }
    
    printf("rate of acceptation : %lf\n", (double) *count / (double) (T_iter*(N_eq + N_tot)));

    char filename1[64];
    sprintf(filename1, "iterations/sqr_avg_mag_vsT_L%d.txt", Lx);
    write_to_file_iter(sqr_avg_mag_vsT, T_iter, filename1);

    char filename2[64];
    sprintf(filename2, "iterations/avg_E_vsT_L%d.txt", Lx);
    write_to_file_iter(avg_E_vsT, T_iter, filename2);


    free(model->lat);
    free(model);
    free(state);
    free(sqr_avg_mag_arr);
    free(sqr_avg_mag_vsT);
    free(avg_E_arr);
    free(avg_E_vsT);
    return 0;
}


