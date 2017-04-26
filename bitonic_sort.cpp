#include <omp.h>
#include <cstdlib>
#include <cstdio>
#include <iostream>

using namespace std;

// constants
#define UP true
#define DOWN false

//pre declarations
void print_array(int *arr, int n);

void read_array(int *arr, int n);


void compareAndSwap(int *seq, int seq_length, bool dir) {

    int middle = seq_length / 2;

    // compare & swap a[i] <-> a[i + n/2]
    #pragma omp parallel for
    for (int i = 0; i < middle; i++) {
        if (dir == (seq[i] > seq[i + middle])) {
            swap(seq[i], seq[i + middle]);
        }
    }
}

// sort bitonic sequence in 'dir' order
void merge(int *seq, int seq_length, bool dir) {

    if (seq_length == 1) return;

    compareAndSwap(seq, seq_length, dir);

    int middle = seq_length / 2;

    #pragma omp parallel for
    for (int i = 0; i < 2; i++) {
        if (i == 0) merge(seq, middle, dir);
        else merge(&seq[middle], seq_length - middle, dir);
    }
}


// FIXME: support odd
void bitonicSort(int *seq, int seq_length, bool dir) {

    if (seq_length == 1) return;

    // recursive calls in different threads
    #pragma omp parallel for
    for (int i = 0; i < 2; i++) {
        if (i == 0) bitonicSort(seq, seq_length / 2, UP); // ASC order
        else bitonicSort(&seq[seq_length / 2], seq_length - seq_length / 2, DOWN); // DESC order
    }

    merge(seq, seq_length, dir);
}


int main() {
    int n;
    cin >> n;

    int *arr = new int[n];
    read_array(arr, n);

    cout << "Entered array:\n";
    print_array(arr, n);

    bitonicSort(arr, n, UP);

    cout << "Sorted array:\n";
    print_array(arr, n);

    delete[] arr;

    return 0;
}


void read_array(int *arr, int n) {
    for (int i = 0; i < n; i++) {
        cin >> arr[i];
    }
}

void print_array(int *arr, int n) {
    cout << "[ ";
    for (int i = 0; i < n; i++) {
        cout << arr[i] << " ";
    }
    cout << "]" << endl;
}
