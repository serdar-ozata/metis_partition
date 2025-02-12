//
// Created by serdar on 1/2/25.
//

#ifndef METIS_PARTITION_TYPEDEF_H
#define METIS_PARTITION_TYPEDEF_H
#include <cstdlib>
#include <parhip_interface.h>

typedef struct {
    double weight;
    int out;
} out_edge;


typedef struct {
    int rows, total_rows;
    idxtype nnz;
    idxtype* xadj;
    idxtype* adjncy;
    idxtype* adjwgt;
    idxtype* vwgt;
    idxtype* vtxdist;
} SparseMat;

void deleteSparseMat(SparseMat &m);

void printMat(SparseMat &m, bool serial);

void printAvailMem();

#endif //METIS_PARTITION_TYPEDEF_H
