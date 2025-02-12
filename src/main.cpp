#include <iostream>
#include "io.h"
#include "partition.h"
#include <mpi.h>
#include <algorithm>
#include "parhip_interface.h"
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
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) printAvailMem();
    SparseMat mat = readFile(input_file);
    idxtype *partition = partitionWithMetis(mat, n_parts);
    writeInpart(input_file, output_path, partition, n_parts, mat);

    // free memory
    delete[] partition;
    MPI_Finalize();

    return 0;
}
