#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <math.h>
#include <stdbool.h>  
#include "../utils.h"
#include "../lattice.h"
#include "../model.h"

int main(){
    SquareLattice* lat = malloc(sizeof(SquareLattice));
    lat->L = 9;
    lat->Lx = 3;
    lat->Ly = 3;
    lat->pbc = false;

    ising* model = define_model(1, 0, 1, lat);
    int* state = malloc(lat->L * sizeof(int));
    for(int i = 0; i < lat->L; i++){
        state[i] = 1;
    }

    int* state2 = malloc(lat->L * sizeof(int));
    for(int i = 0; i < lat->L; i++){
        if (i % 2 == 0){
            state2[i] = 1;
        } else {
            state2[i] = -1;
        }
    }

    printf("total energy for all the spin up : %lf\n", total_energy(model, state));
    printf("total energy for spin up and down : %lf\n", total_energy(model, state2));
}