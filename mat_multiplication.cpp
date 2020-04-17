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
            printf("%f ", a[i * N1 + j]);
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

void init_sizes_of_parts(int *dims, int N1, int N3, int *parts_row_size, int *parts_column_size)
{
    int row_size = N1 / dims[0];              //размер каждой вертикальной полосы
    int column_size = N3 / dims[1];           //размер каждой горизонтальной полосы

    int row_remainder_part = N1 % dims[0];       //размер оставшейся вертикальной полосы
    int column_remained_part = N3 % dims[1];     //размер оставшейся горизонтальной полосы

    for (int i = 0; i < dims[0]; ++i)
    {
        parts_row_size[i] = row_size;
    }

    for (int i = 0; i < dims[1]; ++i)
    {
        parts_column_size[i] = column_size;
    }

    parts_row_size[dims[0] - 1] += row_remainder_part;
    parts_column_size[dims[1] - 1] += column_remained_part;
}

void init_property_of_parts_A(int *dims, int N2, int *parts_row_size, int *parts_A_size, int *parts_A_displs)
{
    for (int i = 0; i < dims[0]; ++i)
    {
        parts_A_size[i] = parts_row_size[i] * N2;
        parts_A_displs[i] = parts_row_size[i] * N2 * i;
    }
}

void init_property_of_parts_B(int *dims, int *parts_B_size, int *parts_B_displs)
{
    for (int i = 0; i < dims[1]; ++i)
    {
        parts_B_displs[i] = i;
        parts_B_size[i] = 1;
    }
}

void mat_multiplication(double *a, double *b, double *c, int N1, int N2, int N3)
{
    int comm_size, comm_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

    int dims[NDIMS] = {0, 0};            //количество процессов по координатам x и y
    MPI_Comm comm_2D;
    create_comm_2D(dims, &comm_2D, comm_size);

    int pr_coords[2];
    MPI_Cart_coords(comm_2D, comm_rank, MAX_DIMS, pr_coords);

    int parts_row_size[dims[0]];        //длина строки для каждого процесса
    int parts_column_size[dims[1]];     //длина столбца для каждого процесса

    init_sizes_of_parts(dims, N1, N3, parts_row_size, parts_column_size);

    ///БУДЕ НУЖНО УДАЛИТЬ ЭТИ ПЕРЕМЕННЫЕ
    int pr_row_size = parts_row_size[0];
    int pr_column_size = parts_column_size[0];

    int row_remainder_part = parts_row_size[dims[0] - 1] - pr_row_size;
    int column_remained_part = parts_column_size[dims[1] - 1] - pr_column_size;

    int current_row_size = parts_row_size[pr_coords[0]];
    int current_column_size = parts_column_size[pr_coords[1]];
    ///ДО ВОТ ЭТОГО МЕСТА

    int *parts_A_displs, *parts_A_size, *parts_B_displs, *parts_B_size;

    MPI_Datatype typeB, typeBMod;

    if (comm_rank == 0)
    {
        parts_A_displs = (int *) calloc(dims[0], sizeof(int));
        parts_A_size = (int *) calloc(dims[0], sizeof(int));
        parts_B_displs = (int *) calloc(dims[1], sizeof(int));
        parts_B_size = (int *) calloc(dims[1], sizeof(int));

        createsTypes(&typeB, &typeBMod, pr_column_size, column_remained_part, N1, N2, N3);

        init_property_of_parts_A(dims, N2, parts_row_size, parts_A_size, parts_A_displs);
        init_property_of_parts_B(dims, parts_B_size, parts_B_displs);
    }

    MPI_Comm comm1DColumns;
    MPI_Comm comm1DRows;
    createComms(comm_2D, &comm1DColumns, &comm1DRows, N1, N2, N3);

    double *aPart = (double *) calloc(current_row_size * N2, sizeof(double));
    double *bPart = (double *) calloc(current_column_size * N2, sizeof(double));
    double *cPart = (double *) calloc(current_column_size * current_row_size, sizeof(double));

    if (pr_coords[1] == 0)
    {
        MPI_Scatterv(a, parts_A_size, parts_A_displs, MPI_DOUBLE, aPart, current_row_size * N2, MPI_DOUBLE, 0,
                     comm1DColumns);
    }

    if (pr_coords[0] == 0)
    {

        MPI_Scatterv(b, parts_B_size, parts_B_displs, typeB, bPart, current_column_size * N2, MPI_DOUBLE, 0,
                     comm1DRows);

        if (pr_coords[1] == 0)
        {
            MPI_Send(b + (N3 / pr_column_size - 1) * pr_column_size, 1, typeBMod, dims[1] - 1, 0, comm1DRows);
        }

        if (dims[1] - 1 == pr_coords[1])
            MPI_Recv(bPart, pr_column_size * N2 + column_remained_part * N2, MPI_DOUBLE, 0, 0, comm1DRows,
                     MPI_STATUS_IGNORE);
    }

    MPI_Bcast(aPart, current_row_size * N2, MPI_DOUBLE, 0, comm1DRows);

    MPI_Bcast(bPart, current_column_size * N2, MPI_DOUBLE, 0, comm1DColumns);

    for (int i = 0; i < current_row_size; ++i)
    {
        for (int j = 0; j < current_column_size; ++j)
        {
            for (int k = 0; k < N2; ++k)
            {
                cPart[i * current_column_size + j] += aPart[i * N2 + k] * bPart[k * current_column_size + j];
            }
        }
    }

    collectCMatrix(dims, pr_coords, c, cPart, pr_row_size, row_remainder_part, pr_column_size, column_remained_part,
                   comm_rank, N1,
                   N2, N3);


    if (comm_rank == 0)
    {
        MPI_Type_free(&typeB);
    }
    free(aPart);
    free(bPart);
    free(cPart);
}

void createComms(MPI_Comm comm2D, MPI_Comm *columns, MPI_Comm *rows, int M, int N, int K)
{
    int remainsRow[2] = {0, 1};
    int remainsColumns[2] = {1, 0};

    MPI_Cart_sub(comm2D, remainsColumns, columns);
    MPI_Cart_sub(comm2D, remainsRow, rows);
}

void createsTypes(MPI_Datatype *typeB, MPI_Datatype *typeBMod, int sizeColumnStrip, int sizeColumnStripMod, int M, int N,
             int K)
{
    //MPI_Type_vector(кол-во блоков, кол-во элементов в блоке, шаг, базовый тип, производный тип);
    MPI_Type_vector(N, sizeColumnStrip, K, MPI_DOUBLE, typeB);
    MPI_Type_vector(N, sizeColumnStrip + sizeColumnStripMod, K, MPI_DOUBLE, typeBMod);

    //MPI_Type_create_resized(*typeB, 0, sizeColumnStrip * sizeof(double), typeB);

    MPI_Type_commit(typeB);
    MPI_Type_commit(typeBMod);

}

void collectCMatrix(int *dims, int *coords, double *c, double *cPart, int sizeRowStrip, int sizeRowStripMod,
                    int sizeColumnStrip, int sizeColumnStripMod, int rank, int M, int N, int K)
{
    int yTmp, xTmp;
    int recvCountsY[dims[0] * dims[1]];
    int recvCountsSplitY[dims[0] * dims[1]];
    int recvCountsX[dims[0] * dims[1]];
    int recvCountsSplitX[dims[0] * dims[1]];
    for (int i = 0; i < dims[0]; i++)
    {
        for (int j = 0; j < dims[1]; j++)
        {
            xTmp = (i == dims[0] - 1) ? (sizeRowStrip + sizeRowStripMod) : sizeRowStrip;
            yTmp = (j == dims[1] - 1) ? (sizeColumnStrip + sizeColumnStripMod) : sizeColumnStrip;
            recvCountsY[i * dims[1] + j] = yTmp * N;
            recvCountsX[i * dims[1] + j] = xTmp * N;
            recvCountsSplitX[i * dims[1] + j] = xTmp;
            recvCountsSplitY[i * dims[1] + j] = yTmp;
        }
    }
    int countA = recvCountsY[coords[1] * dims[1] + coords[1]];
    int splitX = recvCountsSplitX[coords[0] * dims[1] + coords[1]];
    int splitY = recvCountsSplitY[coords[0] * dims[1] + coords[1]];
    sendRecvMatrix(c, cPart, dims, recvCountsY, splitX, splitY, countA, rank, M, N, K);
}

void sendRecvMatrix(double *c, double *cPart, int *dims, int *recvCountsY, int splitX,
                    int splitY, int countA, int rank, int M, int N, int K)
{
    if (rank == 0)
    {
        int index = 0;
        for (int i = 0; i < M * K; i++)
        {

            int trankY = (i % K) / (K / dims[1]);
            //Zero rank gets Y coord for recv and collect (factorization+count the coord)
            int trankX = ((i - (i % K)) / M) / (M / dims[0]);
            //Zero rank gets X coord for recv and collect (factorization+count the coord)

            if (trankX >= dims[0]) trankX = dims[0] - 1;//if the number goes beyond the dims => sub 1
            if (trankY >= dims[1]) trankY = dims[1] - 1;//if the number goes beyond the dims => sub 1

            int trank = trankX * dims[1] + trankY;

            if (trank != 0)
            {

                MPI_Recv(&c[i],
                         recvCountsY[trankX * dims[1] + trankY] / N,
                         MPI_DOUBLE,
                         trankX * dims[1] + trankY,
                         0,
                         MPI_COMM_WORLD,
                         MPI_STATUS_IGNORE);
                i += (recvCountsY[trankX * dims[1] + trankY] / N) - 1;

            } else
            {
                c[i] = cPart[index++];
            }
        }
    } else
    {
        for (int i = 0; i < splitX * splitY; i += countA / N)
        {
            MPI_Send(&cPart[i], splitY, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);//Send a line from subprocess
        }
    }

}
