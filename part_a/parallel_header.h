#include <assert.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* function declaration */
double** init_matrix(int N, int M);
void fill_matrix(double** matrix, int N, int M);
void fill_matrix_with_zero(double** matrix, int N, int M);
double* init_res_arr(int N);
int* init_index_arr(int N);

double inner_product(double* row1, double* row2, int M);
void comp_matrix(double** matrix, int N, int M, double** res);
void comp_arr(double** matrix, int N, int M, double* arr, int* index_arr);

void print_matrix(double** matrix, int N, int M);
void print_upper_tri_matrix(double* arr, int* index_arr, int m_len,
                            int index_len);

/* function initalization */

/* Intialize a matrix with N rows and M columns */
double** init_matrix(int N, int M) {
    double* flatten = (double*)malloc(N * M * sizeof(double));
    double** matrix = (double**)malloc(N * sizeof(double*));
    int i;

    memset(flatten, 0, N * M * sizeof(double));

    for (i = 0; i < N; i++) {
        matrix[i] = &flatten[i * M];
    }

    return matrix;
}

void fill_matrix(double** matrix, int N, int M) {
    int i, j;
    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            matrix[i][j] = ((double)rand() / RAND_MAX * 100.0);
        }
    }
}

void fill_matrix_with_zero(double** matrix, int N, int M) {
    int i, j;
    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            matrix[i][j] = 0;
        }
    }
}

double* init_res_arr(int N) {
    double* arr = (double*)malloc(N * sizeof(double));
    int i;

    for (i = 0; i < N; i++) {
        arr[i] = 0;
    }

    return arr;
}

int* init_index_arr(int N) {
    int* arr = (int*)malloc(N * sizeof(int));
    int i, row;

    for (i = 0, row = N; i < N; i++, row--) {
        arr[i] = (i == 0) ? row : row + arr[i - 1];
    }

    return arr;
}

double inner_product(double* row1, double* row2, int M) {
    float result = 0;
    int i;

    for (i = 0; i < M; i++) {
        result += row1[i] * row2[i];
    }

    return result;
}

void comp_matrix(double** matrix, int N, int M, double** res) {
    int i, j;

    for (i = 0; i < N; i++) {
        for (j = i; j < N; j++) {
            res[i][j - i] = inner_product(matrix[i], matrix[j], M);
        }
    }
}

void comp_arr(double** matrix, int N, int M, double* arr, int* index_arr) {
    int i, j;
    int cnt;
    cnt = 0;

    for (i = 0; i < N; i++) {
        for (j = i; j < N; j++) {
            arr[cnt] = inner_product(matrix[i], matrix[j], M);
            cnt++;
        }
        // Fill the index_arr
        if (i >= 1) {
            index_arr[i] = index_arr[i - 1] + N - i + 1;
        } else {
            index_arr[i] = 0;
        }
    }
    // Fill the index_arr when i == N - 1
    index_arr[i] = index_arr[i - 1] + N - i + 1;
}

void print_matrix(double** matrix, int N, int M) {
    int i, j;
    printf("matrix: \n");
    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            printf("%.2f ", matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void print_upper_tri_matrix(double* arr, int* index_arr, int m_len,
                            int index_len) {
    int i, j;

    for (i = 0, j = 0; i <= index_len; i++) {
        int curr_len;
        curr_len = index_arr[i];
        for (; j < curr_len; j++) {
            printf("%.2f ", arr[j]);
        }
        printf("\n");
    }
}