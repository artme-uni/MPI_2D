#include <iostream>
#include "mat_multiplication.h"

void mat_init(double *a, double *b, int N1, int N2, int N3)
{
    for (int i = 0; i < N1; ++i)
    {
        for (int j = 0; j < N2; ++j)
        {
            a[i * N2 + j] = 1;
        }
    }

    for (int i = 0; i < N2; ++i)
    {
        for (int j = 0; j < N3; ++j)
        {
            b[i * N3 + j] = 1;
        }
    }
}

void mat_print(double *a, int N1, int N2)
{
    for (int i = 0; i < N1; ++i)
    {
        for (int j = 0; j < N2; ++j)
        {
            printf("%f ", a[i * N2 + j]);
        }
        printf("\n");
    }
}

void create_comm_2D(int *dims, MPI_Comm *comm_2D, int comm_size)
{
    int periods[NDIMS] = {FALSE, FALSE};    //периодичность решетки в каждой размерности
    int reorder = FALSE;                    //переупорядочивание нумерации

    MPI_Dims_create(comm_size, NDIMS, dims);

    MPI_Cart_create(MPI_COMM_WORLD, NDIMS, dims, periods, reorder, comm_2D); //конструктор
}

void createTypeB(MPI_Datatype *typeB, int columns_count, int N2, int N3)
{
    //MPI_Type_vector(кол-во блоков, кол-во элементов в блоке, шаг, базовый тип, производный тип);
    MPI_Type_vector(N2, columns_count, N3, MPI_DOUBLE, typeB);
    MPI_Type_create_resized(*typeB, 0, columns_count * sizeof(double), typeB);
    MPI_Type_commit(typeB);
}

void createTypeC(MPI_Datatype *typeC, int columns_count, int rows_count, int N3)
{
    //MPI_Type_vector(кол-во блоков, кол-во элементов в блоке, шаг, базовый тип, производный тип);
    MPI_Type_vector(rows_count, columns_count, N3, MPI_DOUBLE, typeC);
    MPI_Type_create_resized(*typeC, 0, columns_count * sizeof(double), typeC);
    MPI_Type_commit(typeC);
}

void create_sub_comms(MPI_Comm comm_2D, MPI_Comm *comm_columns, MPI_Comm *comm_rows)
{
    int rows[2] = {0, 1};
    int columns[2] = {1, 0};

    MPI_Cart_sub(comm_2D, columns, comm_columns);
    MPI_Cart_sub(comm_2D, rows, comm_rows);
}

void init_parts_A_value(double *parts_A_value, double *A, int parts_A_size, int *pr_coords, MPI_Comm comm_Column,
                        MPI_Comm comm_Row)
{
    if (pr_coords[1] == 0)
    {
        MPI_Scatter(A, parts_A_size, MPI_DOUBLE, parts_A_value, parts_A_size, MPI_DOUBLE, 0, comm_Column);
    }

    MPI_Bcast(parts_A_value, parts_A_size, MPI_DOUBLE, 0, comm_Row);
}

void init_parts_B_value(double *parts_B_value, double *B, int parts_B_size, int *pr_coords, MPI_Comm comm_Column,
                        MPI_Comm comm_Row, MPI_Datatype MPI_typeB)
{
    if (pr_coords[0] == 0)
    {
        MPI_Scatter(B, 1, MPI_typeB, parts_B_value, parts_B_size, MPI_DOUBLE, 0, comm_Row);
    }

    MPI_Bcast(parts_B_value, parts_B_size, MPI_DOUBLE, 0, comm_Column);
}

void mult_parts(double *parts_A_value, double *parts_B_value, double *parts_Result_value,
                int rows_count, int columns_count, int N2)
{
    for (int i = 0; i < rows_count; ++i)
    {
        for (int j = 0; j < columns_count; ++j)
        {
            for (int k = 0; k < N2; ++k)
            {
                parts_Result_value[i * columns_count + j] +=
                        parts_A_value[i * N2 + k] * parts_B_value[k * columns_count + j];
            }
        }
    }
}

void init_parts_C_property(int **parts_C_offset, int **parts_C_size, int *dims, int *pr_coords, int rows_count)
{
    if(pr_coords[0] == 0 && pr_coords[1]==0)
    {
        *parts_C_offset = (int *) calloc(dims[0] * dims[1], sizeof(int));
        *parts_C_size = (int *) calloc(dims[0] * dims[1], sizeof(int));

        for (int i = 0; i < dims[0] * dims[1]; ++i)
        {
            (*parts_C_size)[i] = 1;
        }
        for (int i = 0; i < dims[0]; ++i)
        {
            for (int j = 0; j < dims[1]; ++j)
            {
                (*parts_C_offset)[i * dims[1] + j] = i * dims[1] * rows_count + j;
            }
        }
    }
}

int mat_multiplication(MPI_Comm comm_2D, int *dims, int *pr_coords,
                       double *A, double *B, double *Result, int N1, int N2, int N3)
{
    int rows_count = N1 / dims[0];
    int columns_count = N3 / dims[1];

    MPI_Datatype MPI_typeB;
    createTypeB(&MPI_typeB, columns_count, N2, N3);

    MPI_Comm comm_Column;
    MPI_Comm comm_Row;
    create_sub_comms(comm_2D, &comm_Column, &comm_Row);

    int parts_A_size = rows_count * N2;
    int parts_B_size = columns_count * N2;

    double *parts_A_value = (double *) calloc(parts_A_size, sizeof(double));
    double *parts_B_value = (double *) calloc(parts_B_size, sizeof(double));

    init_parts_A_value(parts_A_value, A, parts_A_size, pr_coords, comm_Column, comm_Row);
    init_parts_B_value(parts_B_value, B, parts_B_size, pr_coords, comm_Column, comm_Row, MPI_typeB);

    double *parts_Result_value = (double *) calloc(columns_count * rows_count, sizeof(double));

    mult_parts(parts_A_value, parts_B_value, parts_Result_value, rows_count, columns_count, N2);

    MPI_Datatype MPI_typeC;
    createTypeC(&MPI_typeC, columns_count, rows_count, N3);

    int *parts_C_offset, *parts_C_size;

    init_parts_C_property(&parts_C_offset, &parts_C_size, dims, pr_coords, rows_count);

    MPI_Gatherv(parts_Result_value, columns_count * rows_count, MPI_DOUBLE, Result, parts_C_size, parts_C_offset,
                MPI_typeC, 0, comm_2D);


    free(parts_A_value);
    free(parts_B_value);
    free(parts_Result_value);

    if(pr_coords[0] == 0 && pr_coords[1] == 0)
    {
        free(parts_C_offset);
        free(parts_C_size);
        MPI_Type_free(&MPI_typeB);
    }

    return 0;
}