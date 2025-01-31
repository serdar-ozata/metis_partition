//
// Created by serdar on 1/28/25.
//
#include "metis.h"
#include <iostream>
#include <cstring>
#include <algorithm>
#include "io.h"
#include "quicksort.h"
int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <input_file> <output_path>" << std::endl;
        return 1;
    }
    string filename = argv[1];
    int n, m;
    idx_t nnz;
    idx_t *row_idx, *col_idx;
    bool symmetric = readMMF(filename, n, m, nnz, row_idx, col_idx);
    idx_t *counts = new idx_t[m];
    memset(counts, 0, sizeof(idx_t) * m);
    idx_t *xadj = new idx_t[m + 1];
    idx_t realNEdges = nnz;
    idx_t *adjncy = new idx_t[realNEdges * 2];
    for (idx_t i = 0; i < realNEdges; ++i) {
        idx_t send_vtx = col_idx[i], recv_vtx = row_idx[i];
        counts[send_vtx]++;
        counts[recv_vtx]++;
    }
    xadj[0] = 0;
    for (idx_t i = 0; i < m; ++i) {
        xadj[i + 1] = xadj[i] + counts[i];
    }
    memset(counts, 0, sizeof(idx_t) * m);
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
    printf("Sorting\n");
    if (not symmetric) {
        // remove duplicates
        for (idx_t i = 0; i < m; ++i) {
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
        for (idx_t i = 0; i < m; ++i) {
            idx_t start = xadj[i];
            idx_t end = xadj[i] + counts[i];
            std::copy(adjncy + start, adjncy + end, adjncy + cur_pos);
            std::copy(adjwgt + start, adjwgt + end, adjwgt + cur_pos);
            cur_pos += counts[i];
            xadj[i] = cur_pos - counts[i];
        }
        xadj[m] = cur_pos;
        realNEdges = cur_pos;
        nnz = realNEdges;
    } else {
        // just do the sorting
        for (idx_t i = 0; i < m; ++i) {
            sort(adjncy + xadj[i], adjncy + xadj[i + 1]);
        }
        nnz = 2 * realNEdges;
    }
    delete[] row_idx;
    delete[] col_idx;
    printf("Writing binary file\n");
    // debug
//    for (idx_t i = 0; i < m; ++i) {
//        for (idx_t j = xadj[i]; j < xadj[i + 1]; ++j) {
//            idx_t  other = adjncy[j];
//            if (i < other) {
//                idx_t ret = binary_search(adjncy + xadj[other], adjncy + xadj[other + 1], i);
//                if (ret == -1) {
//                    cout << "Error: Missing edge" << endl;
//                }
//            } else if (i == other) {
//                cout << "Error: Self loop" << endl;
//            } else break;
//        }
//    }
    // debug end
    string out_filename = argv[2];
    if (out_filename.back() != '/') {
        out_filename += '/';
    }
    out_filename += filename.substr(filename.find_last_of('/') + 1) + ".bin";
    writeBinary(out_filename, n, m, nnz, xadj, adjncy, adjwgt);
    delete[] xadj;
    delete[] adjncy;
    delete[] adjwgt;
    return 0;
}