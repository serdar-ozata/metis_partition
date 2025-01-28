//
// Created by serdar on 1/2/25.
//

#ifndef METIS_PARTITION_IO_H
#define METIS_PARTITION_IO_H
#include <string>
#include <vector>
#include "typedef.h"
using namespace std;

bool readMMF(const string &path, int &n, int &m, idx_t &nnz, idx_t *&row_idx, idx_t *&col_idx);

bool writeBinary(const string &path, int n, int m, idx_t nnz, idx_t* xadj, idx_t* adjncy, idx_t* adjwgt);

void writeInpart(const string input_file, const string output_path, const idx_t *partition, const idx_t nParts, SparseMat& mat);

#endif //METIS_PARTITION_IO_H
