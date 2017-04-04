
#include <iostream>
#include <fstream>
#include <cstring>
#include <math.h>


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
        memcpy(matrix + i * n, v, sizeof(float) * n);
    }

    delete[] v;

    return matrix;
}


int main(int argc, char **argv) {

    char* file = (char *) "/Users/Ramil/ClionProjects/open_mpi/matrix";

    int n = 3000;

    float *m = generateMatrix(n);

    writeMatrix(file, m, n);

    delete[] m;

    return 0;
}
