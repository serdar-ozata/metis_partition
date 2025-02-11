//
// Created by serdar on 1/2/25.
//

#include "io.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
bool readMMF(const string &path, idxtype  &n, idxtype  &m, idxtype &nnz, idxtype  *&row_idx, idxtype  *&col_idx) {
    ifstream file(path);
    if (!file) {
        cout << "Error: Unable to open file" << endl;
        exit(1);
    }
    string line;
    getline(file, line);
    bool symmetric = false;
    if (line.find("symmetric") != string::npos) {
        symmetric = true;
    } else if (line.find("general") == string::npos) {
        cout << "Error: DM-Partition only supports general and symmetric matrices" << endl;
    }
    // skip comments
    while (getline(file, line) && (line[0] == '%' || line.empty())) {}
    istringstream iss(line);
    iss >> n >> m >> nnz;
    row_idx = new idxtype [nnz];
    col_idx = new idxtype [nnz];
    idxtype  self_loops = 0;
    for (idxtype  i = 0; i < nnz; i++) {
        getline(file, line);
        istringstream iss(line);
        iss >> row_idx[i] >> col_idx[i];
        col_idx[i]--;
        row_idx[i]--;
        if (row_idx[i] == col_idx[i]) {
            i--;
            nnz--;
            self_loops++;
        }
    }
    file.close();
    printf("self loops: %d\n", self_loops);
    return symmetric;
}

SparseMat readRBSerial(const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        throw runtime_error("Could not open file: " + filename);
    }

    SparseMat matrix;
    std::string line;

    // Step 1: Read the header
    std::getline(file, line);  // Title line (ignored)
    std::getline(file, line);  // Numerical description
    std::istringstream desc(line);
    int ptrLines, idxLines, valLines, rhsLines;
    desc >> ptrLines >> idxLines >> valLines >> rhsLines;

    // Read matrix dimensions and nonzeros
    std::getline(file, line);
    std::istringstream dims(line);
    dims >> matrix.rows >> matrix.total_rows >> matrix.nnz;

    // Step 2: Read the column pointers
    matrix.xadj = new idxtype[matrix.rows + 1];
    for (int i = 0; i < ptrLines; ++i) {
        std::getline(file, line);
        std::istringstream iss(line);
        int value;
        while (iss >> value) {
            matrix.xadj[i] = value - 1;  // Convert 1-based to 0-based index
        }
    }

    // Step 3: Read the row indices
    matrix.adjncy = new idxtype[matrix.nnz];
    for (int i = 0; i < idxLines; ++i) {
        std::getline(file, line);
        std::istringstream iss(line);
        int value;
        while (iss >> value) {
            matrix.adjncy[i] = value - 1;  // Convert 1-based to 0-based index
        }
    }
    file.close();
    return matrix;
}

void writeInpart(const string input_file, const string output_path, const idxtype *partition, const idxtype nParts, SparseMat& mat) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    // gather partitions
    idxtype * all_partitions = nullptr;
    MPI_Request* requests = nullptr;
    idxtype total_rows = mat.total_rows;
    idxtype rows = mat.rows;
    deleteSparseMat(mat);
    if (rank == 0) {
        all_partitions = new idxtype[total_rows];
        requests = new MPI_Request[size - 1];
        int block_size = total_rows / size;
        for (int i = 1; i < size; i++) {
            int displacement = i * block_size;
            int count = i == size - 1 ? total_rows - displacement : block_size;
            MPI_Irecv(all_partitions + displacement, count, MPI_UNSIGNED_LONG_LONG, i, 0, MPI_COMM_WORLD, requests + i - 1);
        }
        copy(partition, partition + rows, all_partitions);
        MPI_Waitall(size - 1, requests, MPI_STATUSES_IGNORE);
    } else {
        MPI_Send(partition, rows, MPI_UNSIGNED_LONG_LONG, 0, 0, MPI_COMM_WORLD);
        return;
    }
    string input_name = input_file.substr(input_file.find_last_of('/') + 1);
    input_name = input_name.substr(0, input_name.find_first_of('.'));
    string output_file = output_path;
    if (output_path.back() != '/') {
        output_file += '/';
    }
    output_file += input_name + ".inpart." + to_string(nParts);
    FILE *file = fopen(output_file.c_str(), "wb");
    if (file == NULL) {
        perror("Error opening file");
        throw runtime_error("Error opening file");
    }
    fwrite(all_partitions, sizeof(idxtype), total_rows, file);
    fclose(file);
    delete[] all_partitions;
    delete[] requests;
}

bool writeBinary(const string &path, idxtype n, idxtype nnz, idxtype* xadj, idxtype* adjncy, idxtype* adjwgt) {
    FILE* file = fopen(path.c_str(), "wb");
    if (file == NULL) {
        perror("Error opening file");
        return false;
    }
    fwrite(&n, sizeof(idxtype), 1, file);
    fwrite(&nnz, sizeof(idxtype), 1, file);
    bool has_weights = adjwgt != NULL;
    fwrite(&has_weights, sizeof(bool), 1, file);
    fwrite(xadj, sizeof(idxtype), n + 1, file);
    fwrite(adjncy, sizeof(idxtype), nnz, file);
    if (adjwgt != NULL) {
        fwrite(adjwgt, sizeof(idxtype), nnz, file);
    }
    fclose(file);
    return true;
}