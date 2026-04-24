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

int move(ising* model, int* state){ // return if accepted or not
    double E0 = total_energy(model, state);
    int site = rand()%model->lat->L; // a random spin in [0; L-1]
    state[site] = -state[site]; // flip the spin
    double En = total_energy(model, state); // first implementation very slow 
    if (En < E0){
        return 1; // control value --> we accepted the move
    } else {
        double r = (double) rand() / (double) RAND_MAX;
        double prop  = boltzman_probability(En - E0, model->T);
        if (r < prop){
            return 1;
        } else {
            state[site] = -state[site]; 
            return 0; // move rejected
        }
    }
}   