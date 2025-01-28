//
// Created by serdar on 1/9/25.
//
#include "typedef.h"
#include <iostream>
using namespace std;
void deleteSparseMat(SparseMat &m) {
    delete[] m.xadj;
    delete[] m.adjncy;
    if (m.adjwgt != NULL)
        delete[] m.adjwgt;
    delete[] m.vwgt;
    delete[] m.vtxdist;
}

void printMat(SparseMat &m, bool serial) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    for (int i = 0; i < size; ++i) {
        if (rank == i) {
            cout << "rows: " << m.rows << ", cols: " << m.total_rows << ", nnz: " << m.nnz << " rank: " << rank << endl;
            cout << "xadj: ";
            for (int j = 0; j <= m.rows; ++j) {
                cout << m.xadj[j] << " ";
            }
            cout << endl;
            cout << "adjncy: ";
            for (int j = 0; j < m.nnz; ++j) {
                cout << m.adjncy[j] << " ";
            }
            cout << endl;
            if (m.adjwgt != NULL) {
                cout << "adjwgt: ";
                for (int j = 0; j < m.nnz; ++j) {
                    cout << m.adjwgt[j] << " ";
                }
                cout << endl;
            }
            cout << "vwgt: ";
            for (int j = 0; j < m.rows; ++j) {
                cout << m.vwgt[j] << " ";
            }
            cout << endl;
            cout << "vtxdist: ";
            for (int j = 0; j <= size; ++j) {
                cout << m.vtxdist[j] << " ";
            }
            cout << endl;
        }
        if (!serial)
            MPI_Barrier(MPI_COMM_WORLD);
    }
}