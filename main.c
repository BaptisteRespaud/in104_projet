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

    // read params from file "config.txt"

    double T, B;
    int Lx, Ly, N_tot, N_eq;

    FILE *fp = fopen("config.txt", "r");
    fscanf(fp, "T %lf \n", &T);
    printf("Temperature T = %lf \n", T);
    
    fscanf(fp, "B %lf \n", &B);
    printf("External field B = %lf \n", B);
    
    fscanf(fp, "Lx %d \n", &Lx);
    fscanf(fp, "Ly %d \n", &Ly);
    printf("Square lattice Lx = %d Ly = %d \n", Lx, Ly);
    
    fscanf(fp, "N_tot %d \n", &N_tot);
    fscanf(fp, "N_eq %d \n", &N_eq);
    printf("Metropolis N_tot = %d \n", N_tot);
    printf("Metropolis N_eq = %d \n", N_eq);

    fclose(fp);

    SquareLattice* lat = malloc(sizeof(SquareLattice));
    lat->L = Lx*Ly; lat->Lx = Lx; lat->Ly = Ly; lat->pbc = true;
    ising* model = define_model(1, B, T, lat);

    int* state = malloc(lat->L * sizeof(int));
    for(int i = 0; i < lat->L; i++){
        double r = (double) rand() / (double) RAND_MAX;
        if (r < 1/2){
            state[i] = -1;
        } else {
            state[i] = 1; 
        }
    }

    int N_iterations = 10000;

    double* avg_magnetization_array = malloc(N_iterations * sizeof(double));
    double* energy_array = malloc(N_iterations * sizeof(double)); 

    for (int i = 0 ; i < N_iterations; i++){
        // compute the energy of the state and the squared magnetization 
        double m = avg_magnetization(model, state);
        avg_magnetization_array[i] = m * m;
        double e = total_energy(model, state);
        energy_array[i] = e;
        int res = move(model, state);
    }

    printf("average squared magnetization after %d iterations and with T = %f: %lf\n", N_iterations, T, calc_mean(avg_magnetization_array, N_iterations));
    printf("average energy after %d iterations and with T = %f: %lf\n", N_iterations, T, calc_mean(energy_array, N_iterations));

    return 0;
}


