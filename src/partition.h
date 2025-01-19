//
// Created by serdar on 1/2/25.
//

#ifndef METIS_PARTITION_PARTITION_H
#define METIS_PARTITION_PARTITION_H
#include <parmetis.h>
#include <cstdlib>
#include <string>
#include "typedef.h"

void idxToCSR(const idx_t *row_idx, const idx_t *col_idx, bool symmetric, SparseMat &m);

SparseMat readFile(const std::string &filename);

idx_t *partitionWithMetis(SparseMat &m, idx_t nParts);

#endif //METIS_PARTITION_PARTITION_H
