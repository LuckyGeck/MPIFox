#pragma once
#include "matrix.h"

namespace NData {
struct TData {
    TData(size_t fullMatrixSize, size_t n)
        : FullMatrixOrder(fullMatrixSize)
        , N(n)
        , LocalA(n)
        , LocalB(n)
        , LocalC(n)
    {}

    TData(TData&& rhs) = default;

    size_t FullMatrixOrder;
    size_t N;
    NMatrix::TSquareMatrix LocalA;
    NMatrix::TSquareMatrix LocalB;
    NMatrix::TSquareMatrix LocalC;
};
} // namespace NData
