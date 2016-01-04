#pragma once
#include <vector>
#include <mpi.h>


namespace NMatrix {

class TSquareMatrix {
public:
    TSquareMatrix(size_t n)
        : N(n)
        , Data(N * N)
    {}

    TSquareMatrix(TSquareMatrix&& rhs) = default;

    size_t Order() const {
        return N;
    }

    double& At(size_t i, size_t j) {
        return Data[i * N + j];
    }

    double At(size_t i, size_t j) const {
        return Data[i * N + j];
    }

    void MultiplyHere(const TSquareMatrix& a, const TSquareMatrix& b);

    // MPI Helpers
    void InitMPIMatrixType();
    void BroadcastMatrix(int root, MPI_Comm comm);
    void SendrecvReplaceMatrix(int source, int dest, MPI_Comm comm);

private:
    size_t N;
    std::vector<double> Data;
    MPI_Datatype MatrixType;
};

} // namespace NMatrix
