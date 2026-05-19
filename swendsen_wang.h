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

struct uf_elem_s {
    int rank ;
    int elem ;
    struct uf_elem_s* parent ;
};
typedef struct uf_elem_s* uf_elem_t;

typedef uf_elem_t* uf_partition_t ;

uf_partition_t uf_initialize(int size);

void uf_free(uf_partition_t tab, int size);

uf_elem_t uf_find(uf_elem_t x);

void uf_union(uf_elem_t a, uf_elem_t b);

bool connect(uf_partition_t p, ising* model, int* state, int site1, int site2);

void swendsen_wang(ising* model, int* state);

void flip_clusters(uf_partition_t p, ising* model, int* state);

#endif