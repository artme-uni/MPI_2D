#include <cstdlib>
#include <cstdio>
#include <mpi.h>

#include "mat_multiplication.h"

#define RATE_N1 1
#define SIZE_N2 5
#define RATE_N3 3

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    double t1 = MPI_Wtime();

    double *mat_A = NULL;           //размерность N1 на N2
    double *mat_B = NULL;           //размерность N2 на N3
    double *mat_result = NULL;      //размерность N1 на N3

    int comm_size, comm_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

    int dims[NDIMS] = {0, 0};       //количество процессов по координатам x и y соответсвенно
    MPI_Comm comm_2D;
    create_comm_2D(dims, &comm_2D, comm_size);

    int size_N1 = RATE_N1 * dims[0];
    int size_N2 = SIZE_N2;
    int size_N3 = RATE_N3 * dims[1];

    int pr_coords[2];
    MPI_Cart_coords(comm_2D, comm_rank, MAX_DIMS, pr_coords);

    if (pr_coords[0] == 0 && pr_coords[1] == 0)
    {
        printf("Comm_2D: %d * %d\n", dims[0], dims[1]);
        printf("N1 = %d, N2 = %d, N3 = %d\n\n", size_N1, size_N2, size_N3);

        mat_A = (double *) calloc(size_N1 * size_N2, sizeof(double));
        mat_B = (double *) calloc(size_N2 * size_N3, sizeof(double));
        mat_result = (double *) calloc(size_N1 * size_N3, sizeof(double));
        mat_init(mat_A, mat_B, size_N1, size_N2, size_N3);
    }

    mat_multiplication(comm_2D, dims, pr_coords, mat_A, mat_B, mat_result, size_N1, size_N2, size_N3);


    double t2 = MPI_Wtime();

    if (comm_rank == 0)
    {
        printf("Time : %f\n", t2 - t1);
        //mat_print(mat_result, size_N1, size_N3);
        free(mat_A);
        free(mat_B);
        free(mat_result);
    }

    MPI_Finalize();
    return 0;
}