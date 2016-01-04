#include "grid.h"

#include <assert.h>
#include <cmath>
#include <mpi.h>

namespace NGrid {
TGridInfo::TGridInfo() {
    MPI_Comm_size(MPI_COMM_WORLD, &TotalWorkers);
    GridOrder = (int)sqrt((double)TotalWorkers);
    // Check that TotalWorkers is perfect square
    assert(GridOrder * GridOrder == TotalWorkers);

    {
        // Prepare entire grid
        int dimensions[2] = {GridOrder, GridOrder};
        int isPeriodic[2] = {1, 1}; // all dims are periodic
        MPI_Cart_create(MPI_COMM_WORLD, 2, dimensions, isPeriodic, 1, &Comm);
    }
    {
        // Getting our rank and coord in grid
        int coords[2];
        MPI_Comm_rank(Comm, &CurRank);
        MPI_Cart_coords(Comm, CurRank, 2, coords);
        CurRow = coords[0];
        CurCol = coords[1];
    }
    {
        // Create row communicator
        int remain_dims[2];
        remain_dims[0] = 0;
        remain_dims[1] = 1;
        MPI_Cart_sub(Comm, remain_dims, &RowComm);

        // Create col communicator
        remain_dims[0] = 1;
        remain_dims[1] = 0;
        MPI_Cart_sub(Comm, remain_dims, &ColComm);
    }
}
} // namespace NGrid
