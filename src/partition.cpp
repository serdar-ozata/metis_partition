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

using namespace std;
idxtype * partitionWithMetis(SparseMat&m, int nParts) {
    int  noEdgeCut = 0;
    idxtype  *partition = new idxtype [m.rows];
    idxtype  flag = 2;
    if (m.adjwgt != NULL) {
        flag |= 1;
    }
    MPI_Comm comm = MPI_COMM_WORLD;
    double imbalance = 1.05;
    ParHIPPartitionKWay(m.vtxdist, m.xadj, m.adjncy, m.vwgt, m.adjwgt, &nParts, &imbalance, false, 0, FASTSOCIAL, &noEdgeCut, partition, &comm);
    // convert to int*
    return partition;
}

void idxToCSR(const idxtype  *row_idx, const idxtype  *col_idx, bool symmetric, SparseMat &m) {
    idxtype  *counts = new idxtype [m.rows];
    memset(counts, 0, sizeof(idxtype ) * m.rows);
    idxtype  *xadj = new idxtype [m.rows + 1];
    ulong realNEdges = m.nnz;
    idxtype  *adjncy = new idxtype [realNEdges * 2];
    for (idxtype  i = 0; i < realNEdges; ++i) {
        idxtype  send_vtx = col_idx[i], recv_vtx = row_idx[i];
        counts[send_vtx]++;
        counts[recv_vtx]++;
    }
    xadj[0] = 0;
    for (idxtype  i = 0; i < m.rows; ++i) {
        xadj[i + 1] = xadj[i] + counts[i];
    }
    memset(counts, 0, sizeof(idxtype ) * m.rows);
    idxtype  *adjwgt = NULL;
    if (not symmetric) {
        adjwgt = new idxtype [realNEdges * 2];
        memset(adjwgt, 0, sizeof(idxtype ) * realNEdges * 2);
    }
    for (idxtype  i = 0; i < realNEdges; ++i) {
        idxtype  send_vtx = col_idx[i], recv_vtx = row_idx[i];
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
        for (idxtype  i = 0; i < m.rows; ++i) {
            quicksort(adjncy, adjwgt, xadj[i], xadj[i + 1] - 1);
            idxtype  j;
            for (j = xadj[i]; j < xadj[i] + counts[i] - 1; ++j) {
                idxtype  recv_vtx = adjncy[j];
                idxtype  next_recv_vtx = adjncy[j + 1];
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
        idxtype  cur_pos = 0;
        for (idxtype  i = 0; i < m.rows; ++i) {
            idxtype  start = xadj[i];
            idxtype  end = xadj[i] + counts[i];
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
        for (idxtype  i = 0; i < m.rows; ++i) {
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
    fread(&m.total_rows, sizeof(idxtype), 1, file);
    idxtype total_nnz;
    fread(&total_nnz, sizeof(idxtype), 1, file);
    bool has_weights = false;
    fread(&has_weights, sizeof(bool), 1, file);
    long adj_start = (static_cast<long>(m.total_rows) + 3) * sizeof(idxtype) + sizeof(bool);
    long adjwgt_start = adj_start + static_cast<long>(total_nnz) * sizeof(idxtype);

    idxtype  block_size = m.total_rows / size;
    long displacement = rank * block_size;
    idxtype  count = rank == size - 1 ? m.total_rows - displacement : block_size;

    m.vtxdist = new idxtype [size + 1];
    for (idxtype  i = 0; i < size; ++i) {
        m.vtxdist[i] = i * block_size;
    }
    m.vtxdist[size] = m.total_rows;
    long currentPos = ftell(file);
    if (rank == 0) {
        printf("Current position: %ld\n", currentPos);
    }
    fseek(file, displacement * sizeof(idxtype ), SEEK_CUR);
    m.xadj = new idxtype [count + 1];
    for (int i = 0; i <= count; ++i) {
        fread(m.xadj + i, sizeof(idxtype ), 1, file);
    }
    fseek(file, adj_start + m.xadj[0] * sizeof(idxtype ), SEEK_SET);
    m.nnz = m.xadj[count] - m.xadj[0];
    m.adjncy = new idxtype [m.nnz];
    for (idxtype  i = 0; i < m.nnz; ++i) {
        fread(m.adjncy + i, sizeof(idxtype ), 1, file);
    }
    if (has_weights) {
        fseek(file, adjwgt_start + m.xadj[0] * sizeof(idxtype ), SEEK_SET);
        m.adjwgt = new idxtype [m.nnz];
        for (idxtype  i = 0; i < m.nnz; ++i) {
            fread(m.adjwgt + i, sizeof(idxtype ), 1, file);
        }
    } else
        m.adjwgt = NULL;
    fclose(file);
    m.rows = count;
//    sleep(5);
    m.vwgt = new idxtype [m.rows];
    idxtype  first = m.xadj[0];
    for (int i = 0; i <= m.rows; ++i) {
        m.xadj[i] -= first;
    }
    m.xadj[m.rows] = m.nnz;
    for (int i = 0; i < m.rows; ++i) {
        m.vwgt[i] = m.xadj[i + 1] - m.xadj[i];
    }
    return m;
}