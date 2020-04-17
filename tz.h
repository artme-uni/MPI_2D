#ifndef MPI_2D_TZ_H
#define MPI_2D_TZ_H

/*
#include <cstdio>
#include <cstdlib>
#include "mat_multiplication.h"

#define NDIMS 2
#define SIZE_N1 4
#define SIZE_N2 4
#define SIZE_N3 5

#define TRUE 1
#define FALSE 0

int main(int argc, char *argv[])
{
    int comm_size, comm_rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

    int dims[NDIMS] = {0, 0};               //количество процессов по координатам x и y
    int periods[NDIMS] = {FALSE, FALSE};    //периодичность решетки в каждой размерности
    int reorder = FALSE;                    //переупорядочивание нумерации

    MPI_Comm comm_2D;                        //коммутатор новой декартовой топологии
    MPI_Dims_create(comm_size, NDIMS, dims);

    MPI_Cart_create(MPI_COMM_WORLD, NDIMS, dims, periods, reorder, &comm_2D); //конструктор

    double *mat_A = NULL;
    double *mat_B = NULL;
    double *mat_result = NULL;

    if (comm_rank == 0)
    {
        mat_result = (double *) calloc(SIZE_N3 * SIZE_N1, sizeof(double));
        mat_init(mat_A, mat_B);
    }

    mat_mult(mat_A, mat_B, mat_result, dims, comm_rank, comm_2D);

    if (comm_rank == 0)
    {
        mat_print(mat_result);
        resources_free(mat_A, mat_B, mat_result);
        MPI_Comm_free(&comm_2D);
    }

    MPI_Finalize();
    return 0;
}*/

#endif //MPI_2D_TZ_H
