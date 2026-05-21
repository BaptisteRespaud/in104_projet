#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <math.h>
#include <stdbool.h> 

#include "lattice.h"
#include "utils.h"
#include "model.h"

double boltzman_probability(double energy, double T){
    return exp(-energy / T);
}

int move(ising* model, int* state){ 
    int L = model->lat->L;
    int site = rand() % L;

    // local energy before flipping 
    double E_loc_before = local_energy(model, state, site);
    // local energy variations if we flip that spin 
    double dE = -2.0 * E_loc_before;

    // Metropolis algorithm
    if (dE <= 0.0) {
        state[site] = -state[site];
        return 1;
    } else {
        double r = (double) rand() / (double) RAND_MAX;
        if (r < exp(-dE / model->T)) {
            state[site] = -state[site];
            return 1;
        } else {
            return 0;
        }
    }
}
