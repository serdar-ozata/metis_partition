//
// Created by serdar on 1/2/25.
//

#ifndef METIS_PARTITION_TYPEDEF_H
#define METIS_PARTITION_TYPEDEF_H
#include <parmetis.h>
#include <cstdlib>
typedef struct {
    double weight;
    int out;
} out_edge;


typedef struct {
    idx_t rows, total_rows;
    idx_t nnz;
    idx_t* xadj;
    idx_t* adjncy;
    idx_t* adjwgt;
    idx_t* vwgt;
    idx_t* vtxdist;
} SparseMat;

void deleteSparseMat(SparseMat &m);

void printMat(SparseMat &m, bool serial);

#endif //METIS_PARTITION_TYPEDEF_H
