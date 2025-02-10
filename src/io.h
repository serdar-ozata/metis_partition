//
// Created by serdar on 1/2/25.
//

#ifndef METIS_PARTITION_IO_H
#define METIS_PARTITION_IO_H
#include <string>
#include <vector>
#include "typedef.h"
using namespace std;

bool readMMF(const string &path, idxtype  &n, idxtype  &m, idxtype &nnz, idxtype *&row_idx, idxtype *&col_idx);

bool writeBinary(const string &path, idxtype n, idxtype nnz, idxtype* xadj, idxtype* adjncy, idxtype* adjwgt);

void writeInpart(const string input_file, const string output_path, const int *partition, const idxtype nParts, SparseMat& mat);

#endif //METIS_PARTITION_IO_H
