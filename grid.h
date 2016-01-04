#pragma once
#include <mpi.h>

namespace NGrid {

struct TGridInfo {
    TGridInfo();

    int TotalWorkers;
    MPI_Comm Comm;     // Talk to whole grid
    MPI_Comm RowComm;  // Talk to row
    MPI_Comm ColComm;  // Task to col

    int GridOrder;

    int CurRow;
    int CurCol;
    int CurRank;
};

} // namespace NGrid
