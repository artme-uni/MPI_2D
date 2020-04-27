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

int mat_multiplication(MPI_Comm comm_2D, int* dims, int* pr_coords, double *A, double *B, double *Result, int N1, int N2, int N3);
void createsTypes(MPI_Datatype *typeB, MPI_Datatype *typeBMod, int sizeColumnStrip, int sizeColumnStripMod, int M, int N, int K);
void fillDataForEachProc(int *dims, int *coords, int *sizeRows, int *sizeCols, int M, int N, int K);
void init_property_of_parts_B(int *dims, int *sendCountsB, int *displsB, int M, int N, int K);
void sendRecvMatrix(double *c, double *CPart,int *dims,int *recvCountsY,int splitX,int splitY,int countA, int rank, int M, int N, int K);
void collectCMatrix(int *dims, int *coords, double *c,
                    double *cPart, int sizeRowStrip, int sizeRowStripMod, int sizeColumnStrip,
                    int sizeColumnStripMod,int rank, int M, int N, int K);
void create_comm_2D(int *dims, MPI_Comm *comm_2D, int comm_size);

#endif //MPI_2D_MAT_MULTIPLICATION_H
