#ifndef SWWANG_H
#define SWWANG_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <math.h>
#include <stdbool.h> 

#include "lattice.h"
#include "utils.h"
#include "model.h"

/* element of the union_find struct
    it contains : 
        - its rank
        - the value
        - its root
*/
struct uf_elem_s {
    int rank ;
    int elem ;
    struct uf_elem_s* parent ;
};
typedef struct uf_elem_s* uf_elem_t;

/* a partition is an array containing all the elements of the union_find struct */
typedef uf_elem_t* uf_partition_t ;

// initializes a uf_partition of size size
uf_partition_t uf_initialize(int size);

// free function of a partition
void uf_free(uf_partition_t tab, int size);

// returns the parent of the elem x
uf_elem_t uf_find(uf_elem_t x);

// makes the union of the class of a and b 
void uf_union(uf_elem_t a, uf_elem_t b);

/* connect site1 and site2 following the swendsen_wang algorithm 
    - never if S_1 != S_2
    - otherwise with probabily 1 - exp (-2J / T)
*/
bool connect(uf_partition_t p, ising* model, int* state, int site1, int site2);

/* flip spins in state of all clusters contained in the union_find partition p */
void flip_clusters(uf_partition_t p, ising* model, int* state);

/* swendsen_wang algorithm */
void swendsen_wang(ising* model, int* state);

#endif