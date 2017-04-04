#include <stdio.h>
#include <fstream>
#include <iostream>
#include <cstring>
#include <mpi.h>
#include <algorithm>
#include <cmath>
#include <vector>


int read_matrix(const char *file, float *&matrix);

void print_vector(const float *vector, const int n);

void print_matrix(const float *matrix, const int n);

void writeMatrix(char *file, float* matrix, int size) ;

float multiply_vectors(const float *vec1, const float *vec2, const int n) {
    float result = 0;
    for (int i = 0; i < n; ++i) {
        result += vec1[i] * vec2[i];
    }
    return result;
}


void multiply(const float *m, float *v, const int n, const int chunk_size, float *&res) {
    float result[chunk_size];
    for (int i = 0; i < chunk_size; ++i) {
        result[i] = multiply_vectors(&m[i * n], v, n);
    }
    memcpy(res, result, sizeof(float) * chunk_size);
}


void init_pageranks(int n, float *&pageranks) {
    std::fill(pageranks, pageranks + n, 1.f / n);
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

int main(int argc, char **argv) {

    double start, end;

    const int iterations = 1000;

    int rank = 0; // ранк текущего процесса
    int world_size = 0;

    int n = 0; // размер матрицы n x n
    int chunk_size = 0; //
    int left_chunk_size = 0;

    float *matrix = nullptr;
    float *pageranks = nullptr;
    float *sub_matrix = nullptr;
    float *sub_pageranks = nullptr;

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        n = read_matrix(argv[1], matrix);
        left_chunk_size = n % world_size;
    }

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    chunk_size = n / world_size;

    pageranks = new float[n];
    sub_matrix = new float[chunk_size * n];
    sub_pageranks = new float[chunk_size];

    if (rank == 0) {
        init_pageranks(n, pageranks);
    }

    if (rank == 0) {
        MPI_Scatter(matrix, chunk_size * n, MPI_FLOAT, sub_matrix, chunk_size * n, MPI_FLOAT, 0, MPI_COMM_WORLD);
    } else {
        // приняли часть матрицы
        MPI_Scatter(sub_matrix, chunk_size * n, MPI_FLOAT, sub_matrix, chunk_size * n, MPI_FLOAT, 0, MPI_COMM_WORLD);
    }

    for (int i = 0; i < iterations; ++i) {

        MPI_Bcast(pageranks, n, MPI_FLOAT, 0, MPI_COMM_WORLD);

        multiply(sub_matrix, pageranks, n, chunk_size, sub_pageranks);

        if (rank == 0 && left_chunk_size > 0) {

            float left_pagerank[left_chunk_size];

            for (int chunk = left_chunk_size; chunk > 0; --chunk) {

                left_pagerank[left_chunk_size - chunk] = multiply_vectors(&matrix[n * (n - chunk)], pageranks, n);
            }

            for (int j = 0; j < left_chunk_size; ++j) {
                pageranks[n - (left_chunk_size - j)] = left_pagerank[j];
            }
        }

        MPI_Gather(sub_pageranks, chunk_size, MPI_FLOAT, pageranks, chunk_size, MPI_FLOAT, 0, MPI_COMM_WORLD);

        if (rank == 0) {
            normalize_pageranks(pageranks, n);
        }
    }

    if (rank == 0) {
//        print_vector(pageranks, n);
    }

    delete[] pageranks;
    delete[] sub_pageranks;
    delete[] matrix;
    delete[] sub_matrix;


    MPI_Finalize();

    return 0;
}


float *generateMatrix(size_t n) {

    srand((unsigned int) time(0));

    float *matrix = new float[n * n];
    float *v = new float[n];

    memset(v, 0, sizeof(float) * n);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j) {
                v[j] = 0;
                continue;
            }
            v[j] = std::floor((float) rand() / (RAND_MAX / n));
        }
        normalize_pageranks(v, (int) n);
        memcpy(matrix + i * n, v, sizeof(float) * n);
    }

    delete[] v;

    return matrix;
}


void writeMatrix(char *file, float* matrix, int size) {

    std::ofstream fos(file);

    fos << size << "\n";

    for (int i = 0; i < size; ++i) {

        for (int j = 0; j < size; ++j) {

            fos << matrix[i * size + j] << "\t\t";
        }

        fos << std::endl;
    }

    fos.close();
}


void print_matrix(const float *matrix, const int n) {

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            printf("%12.6f", matrix[i * n + j]);
        }
        printf("\n");
    }
}

void print_vector(const float *vector, const int n) {

    printf("( ");
    for (int i = 0; i < n - 1; ++i) {
        printf("%4.6f\t", vector[i]);
    }
    printf("%0.6f )\n", vector[n - 1]);
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



