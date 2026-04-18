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

double ext_field_energy(ising *model, int *state, int site);

double pair_energy(ising *model, int *state, int site1, int site2);

double local_energy(ising *model, int *state, int site);

double total_energy(ising *model, int *state);

double avg_magnetization(ising *model,int *state);



#endif