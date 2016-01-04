#include "matrix.h"
#include <mpi.h>


namespace NMatrix {

void TSquareMatrix::BroadcastMatrix(int root, MPI_Comm comm) {
    MPI_Bcast(this, 1, MatrixType, root, comm);
}

void TSquareMatrix::SendrecvReplaceMatrix(int source, int dest, MPI_Comm comm) {
    MPI_Status status;
    MPI_Sendrecv_replace(this, 1, MatrixType, dest, 0, source, 0, comm, &status);
}

void TSquareMatrix::InitMPIMatrixType() {
    MPI_Datatype  temp_mpi_t;
    MPI_Type_contiguous(N * N, MPI_DOUBLE, &temp_mpi_t);

    int block_lengths[2] = {1, 1};
    MPI_Aint displacements[2];
    {
        MPI_Aint objectStart;
        MPI_Aint address;
        MPI_Address(this, &objectStart);
        MPI_Address(&this->N, &address);
        displacements[0] = address - objectStart;
        MPI_Address(&this->Data[0], &address);
        displacements[1] = address - objectStart;
    }
    MPI_Datatype typelist[2] = {MPI_INT, temp_mpi_t};
    MPI_Type_struct(2, block_lengths, displacements, typelist, &MatrixType);
    MPI_Type_commit(&MatrixType);
}

void TSquareMatrix::MultiplyHere(const TSquareMatrix& a, const TSquareMatrix& b) {
    for (size_t i = 0; i < a.N; ++i) {
        for (size_t j = 0; j < a.N; ++j) {
            for (size_t k = 0; k < b.N; ++k) {
                At(i, j) += a.At(i, k) * b.At(k, j);
            }
        }
    }
}

} // namespace NMatrix
