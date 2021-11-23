#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <stdint.h>
#include <inttypes.h>

void get_chunk(int a, int b, int commsize, int rank, int *lb, int *ub)
{
    int n = b - a + 1;
    int q = n / commsize;

    if (n % commsize)
        q++;

    int r = commsize * q - n;

    int chunk = q;

    if (rank >= commsize - r)
        chunk = q - 1;

    *lb = a;

    if (rank > 0) {
        if (rank <= commsize - r)
            *lb += q * rank;
        else
            *lb += q * (commsize - r) + (q - 1) * (rank - (commsize - r));
    }

    *ub = *lb + chunk - 1;
}

void sgemv(float *a, float *b, float *c, int m, int n, int commsize, int rank, int lb, int ub, int nrows)
{
    for (int i = 0; i < nrows; i++) {
        c[lb + i] = 0.0;
        for (int j = 0; j < n; j++)
            c[lb + i] += a[i * n + j] * b[j];
    }

    if (rank == 0) {
        int *displs = malloc(sizeof(*displs) * commsize);
        int *rcounts = malloc(sizeof(*rcounts) * commsize);

        for (int i = 0; i < commsize; i++) {
            int l, u;
            get_chunk(0, m - 1, commsize, i, &l, &u);
            rcounts[i] = u - l + 1;
            displs[i] = (i > 0) ? displs[i - 1] + rcounts[i - 1]: 0;
        }

        MPI_Gatherv(MPI_IN_PLACE, ub - lb + 1, MPI_FLOAT, c, rcounts, displs, MPI_FLOAT, 0,
                MPI_COMM_WORLD);
    } else {
        MPI_Gatherv(&c[lb], ub - lb + 1, MPI_FLOAT, NULL, NULL, NULL, MPI_FLOAT, 0, MPI_COMM_WORLD);
    }
}

int main(int argc, char **argv)
{
    int rank, commsize;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int n = 45000;
    int m = n;

    double t = MPI_Wtime();

    int lb, ub;
    get_chunk(0, m - 1, commsize, rank, &lb, &ub);
    int nrows = ub - lb + 1;

    float *a = malloc(sizeof(*a) * nrows * n);
    float *b = malloc(sizeof(*b) * n);
    float *c = malloc(sizeof(*c) * m);

    for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < n; j++) {
            a[i * n + j] = lb + i + 1;
        }
    }

    for (int j = 0; j < n; j++) {
        b[j] = j + 1;
    }

    sgemv(a, b, c, m, n, commsize, rank, lb, ub, nrows);

    t = MPI_Wtime() - t;

    double t_max = 0.0;
    MPI_Reduce(&t, &t_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        for (int i = 0; i < m; i++) {
            float r = (i + 1) * (n / 2.0 + pow(n, 2) / 2.0);
            if (fabs(c[i] - r) > 1E-6) {
                fprintf(stderr, "Validation failed: elem %d = %f (real value %f)\n", i, c[i], r);
                break;
            }
        }

        printf("SGEMV: matrix-vector product (c[m] = a[m, n] * b[n]; m = %d, n = %d)\n", m, n);
        printf("Memory used: %" PRIu64 " MiB\n", (uint64_t)(((float)m * n + m + n) * sizeof(double)) >> 20);
        float gflop = 2.0 * m * n * 1E-9;
        printf("Elapsed time (%d procs): %.6f sec.\n", commsize, t_max);
        printf("Performance: %.2f GFLOPS\n", gflop / t_max);

        // printf("%.6f\n", t);
    }

    free(a);
    free(b);
    free(c);

    MPI_Finalize();

    return 0;
}
