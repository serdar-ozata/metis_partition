//
// Created by serdar on 1/2/25.
//

#ifndef METIS_PARTITION_PARTITION_H
#define METIS_PARTITION_PARTITION_H
#include <cstdlib>
#include <string>
#include "typedef.h"
#include "parhip_interface.h"

void idxToCSR(const idxtype *row_idx, const idxtype  *col_idx, bool symmetric, SparseMat &m);

SparseMat readFile(const std::string &filename);

int  *partitionWithMetis(SparseMat &m, int  nParts);

#endif //METIS_PARTITION_PARTITION_H
