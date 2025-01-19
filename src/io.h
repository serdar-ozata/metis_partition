//
// Created by serdar on 1/2/25.
//

#ifndef METIS_PARTITION_IO_H
#define METIS_PARTITION_IO_H
#include <string>
#include <vector>
#include "parmetis.h"
#include "typedef.h"
using namespace std;

enum FILE_TYPE {
    MMF,
    RB
};

bool readMMF(const string &path, int &n, int &m, ulong &nnz, idx_t *&row_idx, idx_t *&col_idx);

void writeInpart(const string input_file, const string output_path, const idx_t *partition, const idx_t nParts, SparseMat& mat);

#endif //METIS_PARTITION_IO_H
