#include "fox.h"
#include <vector>

namespace NFox {
void FoxAlgorithm(NGrid::TGridInfo& grid, NData::TData& data) {
    // Calculate addresses (source - down neighbour; dest - up neighbour)
    int sourceNeighbour = (grid.CurRow + 1) % grid.GridOrder;
    int destNeighbour = (grid.CurRow + grid.GridOrder - 1) % grid.GridOrder;

    NMatrix::TSquareMatrix tempMatrix(data.N);
    tempMatrix.InitMPIMatrixType();

    for (size_t stage = 0; stage < grid.GridOrder; ++stage) {
        int curStageRoot = (grid.CurRow + stage) % grid.GridOrder;
        if (curStageRoot == grid.CurCol) {
            data.LocalA.BroadcastMatrix(curStageRoot, grid.RowComm);
            data.LocalC.MultiplyHere(data.LocalA, data.LocalB);
        } else {
            tempMatrix.BroadcastMatrix(curStageRoot, grid.ColComm);
            data.LocalC.MultiplyHere(tempMatrix, data.LocalB);
        }
        data.LocalB.SendrecvReplaceMatrix(sourceNeighbour, destNeighbour, grid.ColComm);
    }

}
} // namespace NFox
