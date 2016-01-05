#include "fox.h"
#include "grid.h"
#include "matrix.h"
#include "data.h"

#include <mpi.h>

#include <fstream>
#include <iostream>
#include <vector>


namespace {
void ReadMatrixAndScatter(size_t matrixOrder, size_t partMatrixOrder,
                          std::ifstream& input, NGrid::TGridInfo& grid,
                          NMatrix::TSquareMatrix& matrix)
{
    std::vector<double> tmpRow(partMatrixOrder);
    for (int matrixRow = 0; matrixRow < matrixOrder; ++matrixRow) {
        int gridRow = matrixRow / partMatrixOrder;
        int coords[2] = {gridRow, 0};
        for (int gridCol = 0; gridCol < grid.GridOrder; ++gridCol) {
            coords[1] = gridCol;
            int dest = 0;
            MPI_Cart_rank(grid.Comm, coords, &dest);
            if (dest == 0) { // store row into master matrix
                for (int matCol = 0; matCol < partMatrixOrder; ++matCol) {
                    input >> matrix.At(matrixRow, matCol);
                }
            } else { // store tmp row and send it ot dest
                for (auto& elem : tmpRow) {
                    input >> elem;
                }
                MPI_Send(&tmpRow[0], partMatrixOrder, MPI_DOUBLE, dest, 0, grid.Comm);
            }
        }
    }
}

NData::TData ReadInputAndScatter(const char* path, NGrid::TGridInfo& grid) {
    std::ifstream input(path);
    int matrixOrder = 0;
    input >> matrixOrder;
    MPI_Bcast(&matrixOrder, 1, MPI_INT, 0, MPI_COMM_WORLD);
    size_t partMatrixOrder = matrixOrder / grid.GridOrder;
    NData::TData result(matrixOrder, partMatrixOrder);
    ReadMatrixAndScatter(matrixOrder, partMatrixOrder, input, grid, result.LocalA);
    ReadMatrixAndScatter(matrixOrder, partMatrixOrder, input, grid, result.LocalB);
    return result;
}

NData::TData ReceiveMatrix(NGrid::TGridInfo& grid) {
    int matrixOrder = 0;
    MPI_Bcast(&matrixOrder, 1, MPI_INT, 0, MPI_COMM_WORLD);
    size_t partMatrixOrder = matrixOrder / grid.GridOrder;
    NData::TData result(matrixOrder, partMatrixOrder);
    MPI_Status status;
    for (int row = 0; row < partMatrixOrder; ++row) {
        MPI_Recv(&result.LocalA.At(row, 0), partMatrixOrder, MPI_DOUBLE, 0, 0, grid.Comm, &status);
    }
    for (int row = 0; row < partMatrixOrder; ++row) {
        MPI_Recv(&result.LocalB.At(row, 0), partMatrixOrder, MPI_DOUBLE, 0, 0, grid.Comm, &status);
    }
    return result;
}

void DistPrintMatrix(const char* text, NData::TData& data, NMatrix::TSquareMatrix& matrix,
                     NGrid::TGridInfo& grid)
{
    int matrixOrder = data.FullMatrixOrder;
    int partMatrixOrder = data.N;
    if (grid.CurRank == 0) { // i'm a master
        std::cout << text << std::endl;
        std::vector<double> tmpRow(partMatrixOrder);
        for (int matrixRow = 0; matrixRow < matrixOrder; ++matrixRow) {
            int gridRow = matrixRow / partMatrixOrder;
            int coords[2] = {gridRow, 0};
            for (int gridCol = 0; gridCol < grid.GridOrder; ++gridCol) {
                coords[1] = gridCol;
                int source = 0;
                MPI_Cart_rank(grid.Comm, coords, &source);
                if (source == 0) { // load row from master matrix
                    for (int matCol = 0; matCol < partMatrixOrder; ++matCol) {
                        std::cout << matrix.At(matrixRow, matCol) << " ";
                    }
                } else { // store tmp row and send it ot source
                    MPI_Status status;
                    MPI_Recv(&tmpRow[0], partMatrixOrder, MPI_DOUBLE, source, 0, grid.Comm, &status);
                    for (auto& elem : tmpRow) {
                        std::cout << elem << " ";
                    }
                }
            }
            std::cout << std::endl;
        }
    } else { // i'm a slave
        for (int row = 0; row < partMatrixOrder; ++row) {
            MPI_Send(&matrix.At(row, 0), partMatrixOrder, MPI_DOUBLE, 0, 0, grid.Comm);
        }
    }
}

} // anonymous namespace

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int curRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &curRank);
    if (curRank == 0) { // i'm a master
        if (argc < 2) {
            std::cerr << "Usage: mpirun -n 4 " << argv[0] << " <path>" << std::endl;
            return 1;
        }
    }
    NGrid::TGridInfo grid;
    NData::TData data = (curRank == 0) ? ReadInputAndScatter(argv[1], grid)
                                       : ReceiveMatrix(grid);

    DistPrintMatrix("Matrix A:", data, data.LocalA, grid);
    DistPrintMatrix("Matrix B:", data, data.LocalB, grid);

    data.LocalA.InitMPIMatrixType();
    data.LocalB.InitMPIMatrixType();
    NFox::FoxAlgorithm(grid, data);

    DistPrintMatrix("Result:", data, data.LocalC, grid);
    MPI_Finalize();
}
