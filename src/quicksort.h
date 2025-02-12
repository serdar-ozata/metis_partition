//
// Created by serdar on 1/9/25.
//

#ifndef METIS_PARTITION_QUICKSORT_H
#define METIS_PARTITION_QUICKSORT_H
#include <utility>
template <typename T>
void insertion_sort_with_weights(T* adjncy, T* adjwght, T low, T high) {
    for (int i = low + 1; i <= high; ++i) {
        T key = adjncy[i];
        T key_w = adjwght[i];
        int j = i - 1;
        while (j >= low && adjncy[j] > key) {
            adjncy[j + 1] = adjncy[j];
            adjwght[j + 1] = adjwght[j];
            j--;
        }
        adjncy[j + 1] = key;
        adjwght[j + 1] = key_w;
    }
}
template <typename T>
int partition(T* adjncy, T* adjwght, T low, T high) {
    int pivot_idx = (rand() % (high - low + 1)) + low;
    std::swap(adjncy[pivot_idx], adjncy[high]);
    std::swap(adjwght[pivot_idx], adjwght[high]);
    T pivot = adjncy[high];
    int i = low - 1;
    for (int j = low; j < high; ++j) {
        if (adjncy[j] < pivot) {
            i++;
            std::swap(adjncy[i], adjncy[j]);
            std::swap(adjwght[i], adjwght[j]);
        }
    }
    std::swap(adjncy[i + 1], adjncy[high]);
    std::swap(adjwght[i + 1], adjwght[high]);
    return i + 1;
}

template <typename T>
void quicksort(T* adjncy, T* adjwght, T low, T high) {
    const T threshold = 2; // Threshold for small subarrays
    if (low < high) {
        if (high - low + 1 < threshold) {
            insertion_sort_with_weights(adjncy, adjwght, low, high);
        } else {
            T pivotIndex = partition(adjncy, adjwght, low, high);
            quicksort(adjncy, adjwght, low, pivotIndex - 1);
            quicksort(adjncy, adjwght, pivotIndex + 1, high);
        }
    }
}


#endif //METIS_PARTITION_QUICKSORT_H
