//
// Created by serdar on 1/2/25.
//

#include "partition.h"
#include <iostream>
#include <cstring>
#include <algorithm>
#include <mpi.h>
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
    enum FILE_TYPE file_type;
    if (filename.ends_with(".mmf") || filename.ends_with(".mtx")) {
        file_type = MMF;
    } else if (filename.ends_with(".rb")) {
        file_type = RB;
    } else {
        cerr << "Error: Unknown file type" << endl;
        throw runtime_error("Unknown file type");
    }
    switch (file_type) {
        case MMF: {
            SparseMat mat;
            mat.adjwgt = NULL;
            mat.vwgt = NULL;
            mat.xadj = NULL;
            mat.adjncy = NULL;
            int send_counts[size];
            int displs[size];
            memset(send_counts, 0, sizeof(idx_t) * size);
            memset(displs, 0, sizeof(idx_t) * size);
            SparseMat gmat;
            if (rank == 0) {
                int n, m;
                ulong nnz;
                idx_t *row_idx, *col_idx;
                bool symmetric = readMMF(filename, n, m, nnz, row_idx, col_idx);
                gmat.nnz = nnz;
                gmat.rows = n;
                gmat.total_rows = n;
                idxToCSR(row_idx, col_idx, symmetric, gmat);
                idx_t dist_vtx_size = gmat.rows / size;
                for (int i = 0; i < size; ++i) {
                    if (i == size - 1) {
                        send_counts[i] = gmat.rows - i * dist_vtx_size;
                    } else {
                        send_counts[i] = dist_vtx_size;
                    }
                    displs[i] = i * dist_vtx_size;
                }
                mat.total_rows = gmat.total_rows;
//                printMat(gmat, true);
            }
            // vertex weights and pointers
            MPI_Scatter(send_counts, 1, MPI_IDX_T, &mat.rows, 1, MPI_IDX_T, 0, MPI_COMM_WORLD);
            MPI_Bcast(&mat.total_rows, 1, MPI_IDX_T, 0, MPI_COMM_WORLD);
            mat.vwgt = new idx_t[mat.rows];
            MPI_Scatterv(gmat.vwgt, send_counts, displs, MPI_IDX_T, mat.vwgt, mat.rows, MPI_IDX_T, 0, MPI_COMM_WORLD);
            // vtxdist
            mat.vtxdist = new idx_t[size + 1];
            mat.vtxdist[0] = 0;
            if (rank == 0) {
                for (int i = 1; i <= size; ++i) {
                    mat.vtxdist[i] = send_counts[i - 1] + mat.vtxdist[i - 1];
                }
            }
            MPI_Bcast(mat.vtxdist + 1, size, MPI_IDX_T, 0, MPI_COMM_WORLD);

            // xadj
            if (rank == 0) {
                for (int i = 0; i < size; ++i) {
                    send_counts[i]++;
                }
            }
            mat.xadj = new idx_t[mat.rows + 1];
            MPI_Scatterv(gmat.xadj, send_counts, displs, MPI_IDX_T, mat.xadj, mat.rows + 1, MPI_IDX_T, 0, MPI_COMM_WORLD);
            idx_t first_xadj = mat.xadj[0];
            for (int i = 0; i <= mat.rows; ++i) {
                mat.xadj[i] -= first_xadj;
            }
            // edge weights and adjacency
//            if (rank == 0) {
//                for(int i = 0; i < size; ++i) {
//                    send_counts[i] = gmat.xadj[displs[i] + send_counts[i] - 1] - gmat.xadj[displs[i]];
//                }
//                displs[0] = 0;
//                for (int i = 1; i < size; ++i) {
//                    displs[i] = displs[i - 1] + send_counts[i - 1];
//                }
//            }
//            MPI_Scatter(send_counts, 1, MPI_IDX_T, &mat.nnz, 1, MPI_IDX_T, 0, MPI_COMM_WORLD);
//            mat.adjncy = new idx_t[mat.nnz];
//            mat.adjwgt = new idx_t[mat.nnz];
//            MPI_Scatterv(gmat.adjncy, send_counts, displs, MPI_IDX_T, mat.adjncy, mat.nnz, MPI_IDX_T, 0, MPI_COMM_WORLD);
//
//            bool has_adjwgt = gmat.adjwgt != NULL;
//            MPI_Bcast(&has_adjwgt, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
//            if (has_adjwgt) {
//                MPI_Scatterv(gmat.adjwgt, send_counts, displs, MPI_IDX_T, mat.adjwgt, mat.nnz, MPI_IDX_T, 0, MPI_COMM_WORLD);
//            } else {
//                delete[] mat.adjwgt;
//                mat.adjwgt = NULL;
//            }
            mat.nnz = mat.xadj[mat.rows];
            mat.adjncy = new idx_t[mat.nnz];
            bool has_adjwgt = gmat.adjwgt != NULL;
            MPI_Bcast(&has_adjwgt, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
            if (has_adjwgt) {
                mat.adjwgt = new idx_t[mat.nnz];
            } else {
                mat.adjwgt = NULL;
            }
            MPI_Request* request = nullptr;
            #if idx_t == int_32_t
                if (rank == 0) {
                    MPI_Request* send_requests = new MPI_Request[size - 1];
                    for (int i = 1; i < size; ++i) {
                        idx_t* start = gmat.adjncy + gmat.xadj[displs[i]];
                        idx_t* end = gmat.adjncy + gmat.xadj[displs[i] + send_counts[i]];
                        int count = end - start;
                        MPI_Isend(start, count, MPI_IDX_T, i, 0, MPI_COMM_WORLD, send_requests + i - 1);
                    }
                    MPI_Waitall(size, send_requests, MPI_STATUSES_IGNORE);
                    copy(gmat.adjncy, gmat.adjncy + mat.nnz, mat.adjncy);
                } else
                    MPI_Recv(mat.adjncy, mat.nnz, MPI_IDX_T, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            #elif
            // mat.nnz might be larger than int32 max
            int block_count = mat.nnz / INT32_MAX;
            int remaining = mat.nnz % INT32_MAX;

            MPI_Irecv(mat.adjncy + block_count * INT32_MAX, remaining, MPI_IDX_T, 0, 0, MPI_COMM_WORLD, request + block_count);
            if (rank == 0) {
                for (int i = 1; i < size; ++i) {
                    idx_t* start = gmat.adjncy + gmat.xadj[displs[i]];
                    idx_t* end = gmat.adjncy + gmat.xadj[displs[i] + send_counts[i]];
                    unsigned long count = end - start;
                    int block_count = count / INT32_MAX;
                    int remaining = count % INT32_MAX;
                    request = new MPI_Request[block_count + 1];
                    for (int j = 0; j < block_count; ++j) {
                        MPI_Isend(start + j * INT32_MAX, INT32_MAX, MPI_IDX_T, i, 0, MPI_COMM_WORLD, request + j);
                    }
                    MPI_Isend(start + block_count * INT32_MAX, remaining, MPI_IDX_T, i, 0, MPI_COMM_WORLD, request + block_count);
                    MPI_Waitall(block_count + 1, request, MPI_STATUSES_IGNORE);
                    delete[] request;
                }
                copy(gmat.adjncy, gmat.adjncy + mat.nnz, mat.adjncy);
            } else {
                request = new MPI_Request[block_count + 1];
                for (int i = 0; i < block_count; ++i) {
                    MPI_Irecv(mat.adjncy + i * INT32_MAX, INT32_MAX, MPI_IDX_T, 0, 0, MPI_COMM_WORLD, request + i);
                }
                MPI_Waitall(block_count + 1, request, MPI_STATUSES_IGNORE);
                delete[] request;
            }

            #endif
            if (rank == 0) {
                deleteSparseMat(gmat);
            }
            // error check
            if (mat.nnz != mat.xadj[mat.rows]) {
                cerr << "Error: nnz != xadj[rows]" << endl;
                throw runtime_error("Error: nnz != xadj[rows]");
            }
            return mat;
        }
        case RB: {
            throw runtime_error("RB file format is not supported yet");
        }
        default:
            throw runtime_error("Unknown file type");
    }
}