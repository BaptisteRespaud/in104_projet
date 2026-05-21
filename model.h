#ifndef MODEL_H
#define MODEL_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <math.h>
#include <stdbool.h> 
#include "lattice.h"

typedef struct{

    double J;
    double B;
    double T;

    SquareLattice* lat;

}ising;

ising *define_model(double J, double B, double T, SquareLattice *lat);

/* compute the energy from the external magnetic field B according to this formula :
E = - B * spin of the site 
*/
double ext_field_energy(ising *model, int *state, int site);

/* compute the energy from a pair of neighbors according to this formula :
E = - J * S_1 * S_2
*/
double pair_energy(ising *model, int *state, int site1, int site2);

/* compute the local energy of a given site for a state of the system according to this formula :
E = External Field Energy + Each pair energy with the neighbors of the site.
*/
double local_energy(ising *model, int *state, int site);

// compute the total energy of the model in a given state
double total_energy(ising *model, int *state);

/* compute the average magnetization of the system following this formula : 
m = (1 / L) * sum(state[i]) 
*/
double avg_magnetization(ising *model,int *state);



#endif