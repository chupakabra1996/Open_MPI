#include <stdio.h>
#include <fstream>
#include <iostream>
#include <cstring>
#include <mpi.h>
#include <algorithm>


int read_matrix(const char *file, float *&matrix);

void print_vector(const float *vector, int size);


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


void iterate_slave(float *pageranks, int n, int chunk_size, float *sub_pageranks, float *sub_matrix) {

    MPI_Bcast(pageranks, n, MPI_FLOAT, 0, MPI_COMM_WORLD);

    multiply(sub_matrix, pageranks, n, chunk_size, sub_pageranks);

    //NOTE: only root process should have valid buffer
    MPI_Gather(sub_pageranks, chunk_size, MPI_FLOAT, pageranks, chunk_size, MPI_FLOAT, 0, MPI_COMM_WORLD);
}


void iterate_master(float *pageranks, int n, int chunk_size, float *sub_matrix, float *sub_pageranks) {

    MPI_Bcast(pageranks, n, MPI_FLOAT, 0, MPI_COMM_WORLD);

    multiply(sub_matrix, pageranks, n, chunk_size, sub_pageranks);

    MPI_Gather(sub_pageranks, chunk_size, MPI_FLOAT, pageranks, chunk_size, MPI_FLOAT, 0, MPI_COMM_WORLD);

    normalize_pageranks(pageranks, n);
}


void evaluate_slave(int world_size, int iterations) {

    int n; //size of matrix & vector

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD); //get the matrix size

    int chunk_size = n / world_size;

    float *sub_matrix = new float[chunk_size * n];
    float *sub_pageranks = new float[chunk_size];
    float *pageranks = new float[n];

    MPI_Scatter(sub_matrix, chunk_size * n, MPI_FLOAT, sub_matrix, chunk_size * n, MPI_FLOAT, 0, MPI_COMM_WORLD);

    for (int i = 0; i < iterations; ++i) {
        iterate_slave(pageranks, n, chunk_size, sub_pageranks, sub_matrix);
    }

    //clean up
    delete[] sub_matrix;
    delete[] pageranks;
    delete[] sub_pageranks;
}


void evaluate_master(const char *matrix_file, const int world_size, const int iterations) {

    int n; //size of matrix & vector
    float *matrix;

    n = read_matrix(matrix_file, matrix); //read from file

    int chunk_size = n / world_size;

    float *pageranks = new float[n];
    float *sub_matrix = new float[chunk_size * n];
    float *sub_pageranks = new float[chunk_size];

    init_pageranks(n, pageranks); //init it with '1.f/n'

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(matrix, chunk_size * n, MPI_FLOAT, sub_matrix, chunk_size * n, MPI_FLOAT, 0, MPI_COMM_WORLD);

    for (int i = 0; i < iterations; ++i) {
        iterate_master(pageranks, n, chunk_size, sub_matrix, sub_pageranks);
    }
    print_vector(pageranks, n);

    //clean up
    delete[] matrix;
    delete[] sub_matrix;
    delete[] pageranks;
    delete[] sub_pageranks;

}


int main(int argc, char **argv) {

    const int iterations = 1000;
    int rank, world_size;

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) { //master

        evaluate_master(argv[1], world_size, iterations);

    } else { //slaves

        evaluate_slave(world_size, iterations);
    }

    MPI_Finalize();

    return 0;
}


void print_matrix(const float *matrix, const int size) {

    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            printf("%8.3f", matrix[i * size + j]);
        }
        printf("\n");
    }
}

void print_vector(const float *vector, int size) {

    printf("( ");
    for (int i = 0; i < size - 1; ++i) {
        printf("%4.3f\t", vector[i]);
    }
    printf("%0.3f )\n", vector[size - 1]);
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

