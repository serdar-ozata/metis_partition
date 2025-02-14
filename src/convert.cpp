//
// Created by serdar on 1/28/25.
//
#include <iostream>
#include <cstring>
#include <algorithm>
#include "io.h"
#include "quicksort.h"
#include <omp.h>
int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <input_file> <output_path>" << std::endl;
        return 1;
    }
    printAvailMem();
    string filename = argv[1];
    idxtype n, m;
    idxtype nnz;
    idxtype  *row_idx, *col_idx;
    bool symmetric = readMMF(filename, n, m, nnz, row_idx, col_idx);
    printAvailMem();
    idxtype  *counts = new idxtype[m];
    memset(counts, 0, sizeof(idxtype) * m);
    idxtype *xadj = new idxtype [m + 1];
    idxtype realNEdges = nnz;
    idxtype *adjncy = new idxtype [realNEdges * 2];
    for (idxtype i = 0; i < realNEdges; ++i) {
        idxtype send_vtx = col_idx[i], recv_vtx = row_idx[i];
        counts[send_vtx]++;
        counts[recv_vtx]++;
    }
    xadj[0] = 0;
    for (idxtype i = 0; i < m; ++i) {
        xadj[i + 1] = xadj[i] + counts[i];
    }
    memset(counts, 0, sizeof(idxtype) * m);
    idxtype *adjwgt = NULL;
    if (not symmetric) {
        adjwgt = new idxtype [realNEdges * 2];
        memset(adjwgt, 0, sizeof(idxtype) * realNEdges * 2);
    }
    for (idxtype i = 0; i < realNEdges; ++i) {
        idxtype send_vtx = col_idx[i], recv_vtx = row_idx[i];
        adjncy[xadj[send_vtx] + counts[send_vtx]] = recv_vtx;
        adjncy[xadj[recv_vtx] + counts[recv_vtx]] = send_vtx;
        if (not symmetric) {
            adjwgt[xadj[send_vtx] + counts[send_vtx]] = 1;
        }
        counts[send_vtx]++;
        counts[recv_vtx]++;
    }
    cout << "Sorting" << endl;
    if (not symmetric) {
        // remove duplicates
        int local_sorts = 0;
        #pragma omp parallel for firstprivate(local_sorts)
        for (idxtype i = 0; i < m; ++i) {
            quicksort(adjncy, adjwgt, xadj[i], xadj[i + 1] - 1ull);
            idxtype j;
            for (j = xadj[i]; j < xadj[i] + counts[i] - 1; ++j) {
                idxtype recv_vtx = adjncy[j];
                idxtype next_recv_vtx = adjncy[j + 1];
                if (recv_vtx == next_recv_vtx) {
                    counts[i]--;
                    adjwgt[j] |= adjwgt[j + 1]; // if the next one is 1, make it 1
                    // move the remaining elements to one position left
                    size_t remaining = counts[i] + xadj[i] - j - 1;
                    if (remaining > 0) {
                        memmove(adjncy + j + 1, adjncy + j + 2, remaining * sizeof(idxtype));
                        memmove(adjwgt + j + 1, adjwgt + j + 2, remaining * sizeof(idxtype));
                    }
                }
            }
            local_sorts++;
            if (local_sorts % 1000 == 0) {
                cout << "Local sorts: " << local_sorts << endl;
            }
        }
        // now move the unique elements to the beginning
        idxtype cur_pos = 0;
        for (idxtype i = 0; i < m; ++i) {
            idxtype start = xadj[i];
//            std::copy(adjncy + start, adjncy + end, adjncy + cur_pos);
//            std::copy(adjwgt + start, adjwgt + end, adjwgt + cur_pos);
            size_t move_n = counts[i] * sizeof(idxtype);
            memmove(adjncy + cur_pos, adjncy + start, move_n);
            memmove(adjwgt + cur_pos, adjwgt + start, move_n);
            xadj[i] = cur_pos;
            cur_pos += counts[i];
        }
        xadj[m] = cur_pos;
        realNEdges = cur_pos;
        nnz = realNEdges;
    } else {
        // just do the sorting
        #pragma omp parallel for
        for (idxtype i = 0; i < m; ++i) {
            sort(adjncy + xadj[i], adjncy + xadj[i + 1]);
        }
        nnz = 2 * realNEdges;
    }
    delete[] row_idx;
    delete[] col_idx;
    cout << "Writing binary file" << endl;
    // debug
//    for (idxtype i = 0; i < m; ++i) {
//        for (idxtype j = xadj[i]; j < xadj[i + 1]; ++j) {
//            idxtype other = adjncy[j];
//            if (other > ( 1ull << 60)) {
//                cout << "Error: Overflown index" << endl;
//            }
//            if (i < other) {
//                bool ret = binary_search(adjncy + xadj[other], adjncy + xadj[other + 1], i);
//                if (!ret) {
//                    cout << "Error: Missing edge: " << i << " " << other << endl;
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
    writeBinary(out_filename, n, nnz, xadj, adjncy, adjwgt);
    delete[] xadj;
    delete[] adjncy;
    delete[] adjwgt;
    return 0;
}