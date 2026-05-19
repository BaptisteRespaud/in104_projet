#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <math.h>
#include <stdbool.h> 

#include "lattice.h"
#include "utils.h"
#include "model.h"
#include "swendsen_wang.h"

uf_partition_t uf_initialize(int size){
    uf_partition_t tab = malloc(size * sizeof(uf_elem_t));
    for (int i = 0; i < size; i++){
        tab[i] = malloc(sizeof(struct uf_elem_s));
        tab[i]->elem = i;
        tab[i]->rank = 0;
        tab[i]->parent = tab[i];
    }
    return tab;
}

void uf_free(uf_partition_t tab, int size){
    for(int i = 0; i < size; i++){
        free(tab[i]);
    }
    free(tab);
}

uf_elem_t uf_find(uf_elem_t x){
    if (x->parent == x){
        return x;
    } else {
        uf_elem_t representant = uf_find(x->parent);
        x->parent = representant;
        return representant;
    }
}

void uf_union(uf_elem_t a, uf_elem_t b){
    uf_elem_t r_a = uf_find(a);
    uf_elem_t r_b = uf_find(b);
    if (r_a == r_b){
        return;
    } 
    if (r_a->rank < r_b->rank){
        r_a->parent = r_b;
    } else if (r_a->rank > r_b->rank){
        r_b->parent = r_a;
    } else {
        r_a -> parent = r_b;
        r_b->rank += 1;
    }
}

bool connect(uf_partition_t p, ising* model, int* state, int site1, int site2){
    if (state[site1] == state[site2]) {
        double r = (double) rand() / (double) RAND_MAX;
        double prob = 1.0 - exp(-2.0 * model->J / model->T);  
        if (r < prob) {
            uf_union(p[site1], p[site2]);
            return true;
        }
    }
    return false;
}


void flip_clusters(uf_partition_t p, ising* model, int* state){
    int L = model->lat->L;

    int* root_index = malloc(L * sizeof(int));
    for (int i = 0; i < L; i++){
        uf_elem_t root = uf_find(p[i]);
        root_index[i] = root->elem;
    }

    int* cluster_flip = malloc(L * sizeof(int));
    for (int i = 0; i < L; i++){
        cluster_flip[i] = -1;  
    }

    for (int i = 0; i < L; i++){
        int r = root_index[i];
        if (cluster_flip[r] == -1) {
            double u = (double) rand() / (double) RAND_MAX;
            if (u < 0.5){
                cluster_flip[r] = 1;
            } else {
                cluster_flip[r] = 0;
            }
        }
    }

    for (int i = 0 ; i < L; i++){
        int r = root_index[i];
        if (cluster_flip[r]){
            state[i] = -state[i];
        }
    }

    free(root_index);
    free(cluster_flip);
}

void swendsen_wang(ising* model, int* state){

    int Lx = model->lat->Lx;
    int Ly = model->lat->Ly;
    int L  = model->lat->L;
    bool pbc = model->lat->pbc;

    if (!pbc) {
        return;
    }

    uf_partition_t p = uf_initialize(L);

    for (int y = 0; y < Ly; y++){
        for (int x = 0; x < Lx; x++){

            int site = y * Lx + x;

            int right;
            if (x == Lx - 1)
                right = y * Lx;          
            else
                right = site + 1;

            connect(p, model, state, site, right);

            int down;
            if (y == Ly - 1)
                down = x;               
            else
                down = site + Lx;

            connect(p, model, state, site, down);
        }
    }

    flip_clusters(p, model, state);

    uf_free(p, L);
}
