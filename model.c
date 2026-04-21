#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <math.h>
#include <stdbool.h>  
#include "utils.h"
#include "lattice.h"
#include "model.h"


ising* define_model(double J, double B, double T, SquareLattice* lat){

    ising *model = malloc(sizeof(ising));

    model -> J = J;
    model -> B = B;
    model -> T = T;

    model -> lat = lat;

    return model;
}

/* compute the energy from the external magnetic field B according to this formula :
E = - B * spin of the site 
*/
double ext_field_energy(ising *model, int *state, int site){
    return -model->B * state[site];
}

/* compute the energy from a pair of neighbors according to this formula :
E = - J * S_1 * S_2
*/
double pair_energy(ising *model, int *state, int site1, int site2){
    return -model->J * state[site1] * state[site2];
}

/* compute the local energy of a given site for a state of the system according to this formula :
E = External Field Energy + Each pair energy with the neighbors of the site.
*/
double local_energy(ising *model, int *state, int site){
    double ext_field_contribution = ext_field_energy(model, state, site);
    int Lx = model->lat->Lx;
    int Ly = model->lat->Ly;

    if (model->lat->pbc){
        // LEFT
        int left = (site % Lx == 0) ? site + (Lx - 1) : site - 1;
        double left_energy = pair_energy(model, state, site, left);

        // RIGHT
        int right = (site % Lx == Lx - 1) ? site - (Lx - 1) : site + 1;
        double right_energy = pair_energy(model, state, site, right);

        // UP
        int up = (site < Lx) ? site + (Ly - 1) * Lx : site - Lx;
        double up_energy = pair_energy(model, state, site, up);

        // DOWN  (FIXED)
        int down = (site >= (Ly - 1) * Lx) ? site - (Ly - 1) * Lx : site + Lx;
        double down_energy = pair_energy(model, state, site, down);

        return left_energy + right_energy + up_energy + down_energy + ext_field_contribution;
    }
    int row = (int) (site / Lx);
    int col = (int) (site % Lx);

    double pairs_energy = 0;

     if (col > 0)
        pairs_energy += pair_energy(model, state, site, site - 1);

    if (col < Lx - 1)
        pairs_energy += pair_energy(model, state, site, site + 1);

    if (row > 0)
        pairs_energy += pair_energy(model, state, site, site - Lx);

    if (row < Ly - 1)
        pairs_energy += pair_energy(model, state, site, site + Lx);

    return pairs_energy + ext_field_contribution;
}

double total_energy(ising *model, int *state){
    return;   
}