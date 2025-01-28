//
// Created by serdar on 1/2/25.
//

#include "partition.h"
#include <iostream>
#include <cstring>
#include <algorithm>
#include <mpi.h>
#include <csignal>
#include "io.h"
#include "quicksort.h"

#if idx_t == int_32_t
#define MPI_IDX_T MPI_INT
#else
#define MPI_IDX_T MPI_LONG
#endif
using namespace std;
idx_t* partitionWithMetis(SparseMat&m, idx_t nParts) {
    idx_t nCon = 1;  // Number of balancing constraints
    idx_t noEdgeCut = 0;
    idx_t *partition = new idx_t[m.rows];
    int numflag = 0;
    idx_t flag = 2;
    if (m.adjwgt != NULL) {
        flag |= 1;
    }
    MPI_Comm comm = MPI_COMM_WORLD;
    real_t* tpwgts = new real_t[nParts * nCon];
    for (int i = 0; i < nParts; ++i) {
        tpwgts[i] = 1.0 / nParts;
    }
    real_t ubvec = {1.05};
    idx_t options[3] = {0, 0, 0};
    int ret = ParMETIS_V3_PartKway(m.vtxdist, m.xadj, m.adjncy, m.vwgt, m.adjwgt, &flag, &numflag, &nCon, &nParts, tpwgts, &ubvec, options, &noEdgeCut, partition, &comm);

    if (ret != METIS_OK) {
        cout << "METIS failed with error code " << ret << endl;
        exit(1);
    }
    delete[] tpwgts;
    return partition;
}

void idxToCSR(const idx_t *row_idx, const idx_t *col_idx, bool symmetric, SparseMat &m) {
    idx_t *counts = new idx_t[m.rows];
    memset(counts, 0, sizeof(idx_t) * m.rows);
    idx_t *xadj = new idx_t[m.rows + 1];
    ulong realNEdges = m.nnz;
    idx_t *adjncy = new idx_t[realNEdges * 2];
    for (idx_t i = 0; i < realNEdges; ++i) {
        idx_t send_vtx = col_idx[i], recv_vtx = row_idx[i];
        counts[send_vtx]++;
        counts[recv_vtx]++;
    }
    xadj[0] = 0;
    for (idx_t i = 0; i < m.rows; ++i) {
        xadj[i + 1] = xadj[i] + counts[i];
    }
    memset(counts, 0, sizeof(idx_t) * m.rows);
    idx_t *adjwgt = NULL;
    if (not symmetric) {
        adjwgt = new idx_t[realNEdges * 2];
        memset(adjwgt, 0, sizeof(idx_t) * realNEdges * 2);
    }
    for (idx_t i = 0; i < realNEdges; ++i) {
        idx_t send_vtx = col_idx[i], recv_vtx = row_idx[i];
        adjncy[xadj[send_vtx] + counts[send_vtx]] = recv_vtx;
        adjncy[xadj[recv_vtx] + counts[recv_vtx]] = send_vtx;
        if (not symmetric) {
            adjwgt[xadj[send_vtx] + counts[send_vtx]] = 1;
        }
        counts[send_vtx]++;
        counts[recv_vtx]++;
    }
    if (not symmetric) {
        // remove duplicates
        for (idx_t i = 0; i < m.rows; ++i) {
            quicksort(adjncy, adjwgt, xadj[i], xadj[i + 1] - 1);
            idx_t j;
            for (j = xadj[i]; j < xadj[i] + counts[i] - 1; ++j) {
                idx_t recv_vtx = adjncy[j];
                idx_t next_recv_vtx = adjncy[j + 1];
                if (recv_vtx == next_recv_vtx) {
                    counts[i]--;
                    if (adjwgt[j] && adjwgt[j + 1]) {
                        cout << "Error: Duplicate edges with weights" << endl;
                        cout << "send_vtx: " << i << ", recv_vtx: " << recv_vtx <<  ", next_recv_vtx: " << next_recv_vtx << endl;
                    }
                    adjwgt[j] |= adjwgt[j + 1]; // if the next one is 1, make it 1
                    // move the remaining elements to one position left
                    if (j + 2 < xadj[i] + counts[i]) {
                        std::copy(adjncy + j + 2, adjncy + xadj[i + 1], adjncy + j + 1);
                        std::copy(adjwgt + j + 2, adjwgt + xadj[i + 1], adjwgt + j + 1);
                    }

                }
            }
        }
        // now move the unique elements to the beginning
        idx_t cur_pos = 0;
        for (idx_t i = 0; i < m.rows; ++i) {
            idx_t start = xadj[i];
            idx_t end = xadj[i] + counts[i];
            std::copy(adjncy + start, adjncy + end, adjncy + cur_pos);
            std::copy(adjwgt + start, adjwgt + end, adjwgt + cur_pos);
            cur_pos += counts[i];
            xadj[i] = cur_pos - counts[i];
        }
        xadj[m.rows] = cur_pos;
        realNEdges = cur_pos;
        m.nnz = realNEdges;
        int numberOfWeights = 0;
        for (int i = 0; i < realNEdges; ++i) {
            if (adjwgt[i] > 0) {
                numberOfWeights++;
            }
        }
        printf("Number of weights: %d, totalAdj: %d\n", numberOfWeights, realNEdges);
    } else {
        // just do the sorting
        for (idx_t i = 0; i < m.rows; ++i) {
            sort(adjncy + xadj[i], adjncy + xadj[i + 1]);
        }
        m.nnz = 2 * realNEdges;
    }
    m.xadj = xadj;
    m.adjncy = adjncy;
    m.adjwgt = adjwgt;
    m.vwgt = counts;
    delete[] row_idx;
    delete[] col_idx;

}

SparseMat readFile(const std::string &filename) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    SparseMat m;
    FILE *file = fopen(filename.c_str(), "rb");
    if (file == NULL) {
        perror("Error opening file");
        throw runtime_error("Error opening file");
    }
    fread(&m.total_rows, sizeof(int), 1, file);
    fread(&m.rows, sizeof(int), 1, file); // this value is not important
    idx_t total_nnz;
    fread(&total_nnz, sizeof(idx_t), 1, file);
    bool has_weights = false;
    fread(&has_weights, sizeof(bool), 1, file);
    idx_t adj_start = (m.total_rows + 2) * sizeof(idx_t) + 2 * sizeof(int) + sizeof(bool);
    idx_t adjwgt_start = adj_start + total_nnz * sizeof(idx_t);

    int block_size = m.total_rows / size;
    long displacement = rank * block_size;
    int count = rank == size - 1 ? m.total_rows - displacement : block_size;

    m.vtxdist = new idx_t[size + 1];
    for (int i = 0; i < size; ++i) {
        m.vtxdist[i] = i * block_size;
    }
    m.vtxdist[size] = m.total_rows;

    fseek(file, displacement * sizeof(idx_t), SEEK_CUR);
    m.xadj = new idx_t[count + 1];
    for (int i = 0; i <= count; ++i) {
        fread(m.xadj + i, sizeof(idx_t), 1, file);
    }
    fseek(file, adj_start + m.xadj[0] * sizeof(idx_t), SEEK_SET);
    m.nnz = m.xadj[count] - m.xadj[0];
    m.adjncy = new idx_t[m.nnz];
    for (int i = 0; i < m.nnz; ++i) {
        fread(m.adjncy + i, sizeof(idx_t), 1, file);
    }
    if (has_weights) {
        fseek(file, adjwgt_start + m.xadj[0] * sizeof(idx_t), SEEK_SET);
        m.adjwgt = new idx_t[m.nnz];
        for (int i = 0; i < m.nnz; ++i) {
            fread(m.adjwgt + i, sizeof(idx_t), 1, file);
        }
    } else
        m.adjwgt = NULL;
    fclose(file);
    m.rows = count;
//    sleep(5);
    m.vwgt = new idx_t[count];
    idx_t first = m.xadj[0];
    for (int i = 0; i <= count; ++i) {
        m.xadj[i] -= first;
    }
    m.xadj[count] = m.nnz;
    for (int i = 0; i < count; ++i) {
        m.vwgt[i] = m.xadj[i + 1] - m.xadj[i];
    }
    return m;
}