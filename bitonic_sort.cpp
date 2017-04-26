#include <omp.h>
#include <cstdlib>
#include <cstdio>
#include <iostream>

// See http://www.iti.fh-flensburg.de/lang/algorithmen/sortieren/bitonic/oddn.htm

using namespace std;

// Directions
#define ASC true
#define DESC false


//pre declarations
void print_array(int *arr, int n);

void read_array(int *arr, int n);


class BitonicSort {
public:

    // constructor
    BitonicSort(int *seq, int length) :
            seq(seq), length(length) {}

    BitonicSort() : BitonicSort(0, 0) {}

    // sort sequence in some order
    void sort(bool direction) {
        bitonic_sort(0, length, direction);
    }


    void sort() {
        sort(ASC);
    }


    int *get_sequence() const {
        return seq;
    }

    int get_length() const {
        return length;
    }

    // destructor
    ~BitonicSort() {
        if (seq) delete[] seq;
        length = 0;
        seq = 0;
    }

private:
    int *seq; // sequence, array
    int length;


    void bitonic_sort(int low, int seq_length, bool direction) {

        if (seq_length > 1) {

            int middle = seq_length / 2;

            // recursive calls in parallel
            // split sequence into bitonic ones (every sequence of 2 elements is bitonic)
            #pragma omp parallel for
            for (int i = 0; i < 2; ++i) {
                if (i == 0) bitonic_sort(low, middle, !direction);
                else bitonic_sort(low + middle, seq_length - middle, direction);
            }

            // merge bitonic sequences, in other words, sort them considering direction
            merge(low, seq_length, direction);
        }
    }


    // sorts sequence in the 'direction' order
    void merge(int low, int seq_length, bool direction) {

        if (seq_length > 1) {

            // calculate the nearest to the 'seq_length' value, so
            // it's power of 2 ( 2^k = max & max < seq_length )
            int max = max_power_of_two_value_below(seq_length);

            // compare and swap in different threads
            // instead of comparing a[i] with a[i + n/2],
            // we're using a little bit other approach now
            // See the link above
            #pragma omp parallel for
            for (int i = low; i < low + seq_length - max; i++) {
                compare_and_swap(i, i + max, direction);
            }

            // recursive merging in parallel
            #pragma omp parallel for
            for (int i = 0; i < 2; ++i) {
                if (i == 0) merge(low, max, direction);
                else merge(low + max, seq_length - max, direction);
            }

        }
    }


    // compare two sequnce values and swap them if needed
    void compare_and_swap(int i, int j, bool direction) {

        if (direction == (seq[i] > seq[j])) {
            swap(seq[i], seq[j]);
        }
    }

    int max_power_of_two_value_below(int n) {
        int k = 1;
        while (k > 0 && k < n) {
            k <<= 1;
        }
        return k >> 1;
    }
};


/*
 * Usage, i.e `g++ -o bitonic_sort -fopenmp bitonic_sort.cpp`
 *
 * ./bitonic_sort < input.txt
 *
 * where input.txt:
 * 5
 * 5 4 3 2 1
 *
*/
int main() {

    int n;
    cin >> n;

    int *arr = new int[n];
    read_array(arr, n);

    cout << "*** input array ***" << endl;
    print_array(arr, n);

    BitonicSort bitonicSort(arr, n);

    bitonicSort.sort(ASC);

    cout << "*** sorted array ***" << endl;
    print_array(arr, n);

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
