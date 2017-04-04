

#include <cstring>
#include <cstdio>
#include <fstream>

void print_matrix(const float *matrix, const int n) {

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            printf("%12.6f", matrix[i * n + j]);
        }
        printf("\n");
    }
}


//read matrix and return it's size
int read_matrix(const char *file, float *&matrix) {

    std::ifstream fin(file); //create input stream

    int size = 0;
    fin >> size; //read the size of the matrix

    matrix = new float[size * size]; // 2-d matrix in presentation of 1-d vector
    memset(matrix, 0, sizeof(float) * size);

    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            fin >> matrix[j * size + i]; //read matrix values
        }
    }
    fin.close();
    return size;
}

void print_vector(const float *vector, const int n) {

    printf("( ");
    for (int i = 0; i < n - 1; ++i) {
        printf("%4.6f\t", vector[i]);
    }
    printf("%0.6f )\n", vector[n - 1]);
}


void multiply_matrix_to_vector(float *matrix, float*& vector, int n) {

    float res[n]; // result pageranks

    for (int k = 0; k < n; ++k) {
        res[k] = 0.f;
    }

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            res[i] += matrix[i * n + j] * vector[j];
        }
    }

    memcpy(vector, res, sizeof(float) * n);
}


void normalize_pageranks(float *&pageranks, int n) {
    float sum = 0;
    for (int i = 0; i < n; ++i) {
        sum += pageranks[i];
    }
    for (int j = 0; j < n; ++j) {
        pageranks[j] /= sum;
    }
}


void init_pageranks(int n, float *&pageranks) {
    std::fill(pageranks, pageranks + n, 1.f / n);
}


int main(int argc, char** argv) {

    float* matrix = nullptr;
    float* vector = nullptr;

    int size = read_matrix("/Users/Ramil/ClionProjects/open_mpi/matrix",  matrix);
    vector = new float[size];

    init_pageranks(size, vector);

    for (int i = 0; i < 1000; ++i) {

        multiply_matrix_to_vector(matrix, vector, size);

        normalize_pageranks(vector, size);
    }

    delete[] vector;
    delete[] matrix;

    return 0;
}