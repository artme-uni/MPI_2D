#ifndef MPI_2D_MAT_MULTIPLICATION_H
#define MPI_2D_MAT_MULTIPLICATION_H

#include <cstdlib>
#include <cstdio>
#include <mpi.h>

#define MAX_DIMS 2
#define NDIMS 2

#define TRUE 1
#define FALSE 0

void mat_init(double *a, double *b, int N1, int N2, int N3);

void mat_print(double *a, int N1, int N2);

void create_comm_2D(int *dims, MPI_Comm *comm_2D, int comm_size);

void createTypeB(MPI_Datatype *typeB, int columns_count, int N2, int N3);

void createTypeC(MPI_Datatype *typeC, int columns_count, int rows_count, int N3);

void create_sub_comms(MPI_Comm comm_2D, MPI_Comm *comm_columns, MPI_Comm *comm_rows);

void init_parts_A_value(double *parts_A_value, double *A, int parts_A_size, int *pr_coords, MPI_Comm comm_Column,
                        MPI_Comm comm_Row);

void init_parts_B_value(double *parts_B_value, double *B, int parts_B_size, int *pr_coords, MPI_Comm comm_Column,
                        MPI_Comm comm_Row, MPI_Datatype MPI_typeB);

void mult_parts(double *parts_A_value, double *parts_B_value, double *parts_Result_value,
                int rows_count, int columns_count, int N2);
void init_parts_C_property(int **parts_C_offset, int **parts_C_size, int *dims, int *pr_coords, int rows_count);

int mat_multiplication(MPI_Comm comm_2D, int *dims, int *pr_coords,
                       double *A, double *B, double *Result, int N1, int N2, int N3);

#endif //MPI_2D_MAT_MULTIPLICATION_H
