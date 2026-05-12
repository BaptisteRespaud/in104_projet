#ifndef METROPOLIS_H
#define METROPOLIS_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <math.h>
#include <stdbool.h> 

#include "lattice.h"
#include "utils.h"
#include "model.h"


double boltzman_probability(double energy, double T);

int move(ising* model, int* state);   // return if accepted or not



#endif