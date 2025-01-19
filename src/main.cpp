#include <iostream>
#include "io.h"
#include "partition.h"
#include "quicksort.h"
#include <mpi.h>
#include <cstring>
#include <unistd.h>
using namespace std;
int main(int argc, char **argv) {
    if (argc != 4) {
        cerr << "Usage: " << argv[0] << " <input_file> <n_parts> <output_path>" << endl;
        return 1;
    }
    MPI_Init(&argc, &argv);
    string input_file = argv[1];
    int n_parts = stoi(argv[2]);
    string output_path = argv[3];
    SparseMat mat = readFile(input_file);
    //debug
    for (int i = 1; i <= mat.rows; ++i) {
        // check if there is invalid vertex or int overflow
        if (mat.xadj[i] < 0 || mat.xadj[i] > mat.nnz || mat.xadj[i] < mat.xadj[i - 1]) {
            cerr << "Error: invalid vertex or int overflow" << endl;
            throw runtime_error("Error: invalid vertex or int overflow");
        }
        for (int j = mat.xadj[i - 1]; j < mat.xadj[i] - 1; ++j) {
            if (mat.adjncy[j] == mat.adjncy[j + 1]) {
                cerr << "Error: duplicate edge" << endl;
                throw runtime_error("Error: duplicate edge");
            }
        }
    }
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    for (int i = 0; i < size; ++i) {
        if (rank == i) {
            cout << "rows: " << mat.rows << ", cols: " << mat.total_rows << ", nnz: " << mat.nnz << " rank: " << rank << endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    // end debug
    idx_t *partition = partitionWithMetis(mat, n_parts);
    writeInpart(input_file, output_path, partition, n_parts, mat);

    // free memory
    deleteSparseMat(mat);
    delete[] partition;
    MPI_Finalize();

    return 0;
}
