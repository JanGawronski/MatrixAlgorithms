#include "Strassen.h"
#include "SupportFunctions.h"

#include <algorithm>
#include <stdexcept>
#include <memory>
#include <vector>

namespace {

using Matrix = ::Matrix;

Matrix trim(const Matrix &A, int newRows, int newCols) {
    int Arows = static_cast<int>(A.size());
    int Acols = Arows ? static_cast<int>(A[0].size()) : 0;
    if (newRows > Arows || newCols > Acols) {
        throw std::runtime_error("New size must be smaller than original");
    }
    Matrix B = zeroMatrix(newRows, newCols);
    for (int i = 0; i < newRows; i++)
        for (int j = 0; j < newCols; j++)
            B[i][j] = A[i][j];
    return B;
}

Matrix subMatrixPadded(const Matrix &A, int row, int col, int rows, int cols, int paddedRows, int paddedCols) {
    Matrix R = zeroMatrix(paddedRows, paddedCols);
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            R[i][j] = A[row + i][col + j];
    return R;
}

Matrix combine(const Matrix &A11, const Matrix &A12,
               const Matrix &A21, const Matrix &A22) {
    int A11rows = static_cast<int>(A11.size());
    int A12rows = static_cast<int>(A12.size());
    int A11cols = A11rows ? static_cast<int>(A11[0].size()) : 0;
    int A12cols = A12.size() ? static_cast<int>(A12[0].size()) : 0;
    int A21rows = A21.size();
    int A22rows = A22.size();
    Matrix A = zeroMatrix(A11rows + A21rows, A11cols + A12cols);
    for (int i = 0; i < A11rows; i++)
        for (int j = 0; j < A11cols; j++)
            A[i][j] = A11[i][j];

    for (int i = 0; i < A12rows; i++)
        for (int j = 0; j < A12cols; j++)
            A[i][j + A11cols] = A12[i][j];

    for (int i = 0; i < A21rows; i++)
        for (int j = 0; j < A11cols; j++)
            A[i + A11rows][j] = A21[i][j];

    for (int i = 0; i < A22rows; i++)
        for (int j = 0; j < A12cols; j++)
            A[i + A11rows][j + A11cols] = A22[i][j];

    return A;
}

Matrix multiplyRec(const Matrix &A, const Matrix &B) {
    int Arows = static_cast<int>(A.size());
    int Acols = Arows ? static_cast<int>(A[0].size()) : 0;
    int Brows = static_cast<int>(B.size());
    int Bcols = Brows ? static_cast<int>(B[0].size()) : 0;
    if (Arows != Acols || Acols != Brows || Brows != Bcols) {
        throw std::runtime_error("Implemented only for square matrices");
    }

    memCounterEnterCall(static_cast<std::size_t>(Arows), static_cast<std::size_t>(Bcols));

    if (Arows == 1) {
        Matrix Z = A * B;
        memCounterExitCall(static_cast<std::size_t>(Arows), static_cast<std::size_t>(Bcols));
        return Z;
    }

    int A11width = (Acols + 1) / 2;
    int A12width = Acols - A11width;
    int A11height = (Arows + 1) / 2;
    int A21height = Arows - A11height;

    Matrix A11 = subMatrix(A, 0, 0, A11height, A11width);
    Matrix A12 = subMatrixPadded(A, 0, A11width, A11height, A12width, A11height, A11width);
    Matrix A21 = subMatrixPadded(A, A11height, 0, A21height, A11width, A11height, A11width);
    Matrix A22 = subMatrixPadded(A, A11height, A11width, A21height, A12width, A11height, A11width);
    
    Matrix B11 = subMatrix(B, 0, 0, A11height, A11width);
    Matrix B12 = subMatrixPadded(B, 0, A11width, A11height, A12width, A11height, A11width);
    Matrix B21 = subMatrixPadded(B, A11height, 0, A21height, A11width, A11height, A11width);
    Matrix B22 = subMatrixPadded(B, A11height, A11width, A21height, A12width, A11height, A11width);

    Matrix P1 = multiplyRec(A11 + A22, B11 + B22);
    Matrix P2 = multiplyRec(A21 + A22, B11);
    Matrix P3 = multiplyRec(A11, B12 - B22);
    Matrix P4 = multiplyRec(A22, B21 - B11);
    Matrix P5 = multiplyRec(A11 + A12, B22);
    Matrix P6 = multiplyRec(A21 - A11, B11 + B12);
    Matrix P7 = multiplyRec(A12 - A22, B21 + B22);

    Matrix C11 = P1 + P4 - P5 + P7;
    Matrix C12 = trim(P3 + P5, A11height, A12width);
    Matrix C21 = trim(P2 + P4, A21height, A11width);
    Matrix C22 = trim(P1 + P3 - P2 + P6, A21height, A12width);

    Matrix M = combine(C11, C12, C21, C22);

    memCounterExitCall(static_cast<std::size_t>(Arows), static_cast<std::size_t>(Bcols));
    return M;
}

}

class StrassenImpl : public IMnozenie {
public:
    Matrix multiply(const Matrix &A, const Matrix &B) override {
        return multiplyRec(A, B);
    }
};

std::unique_ptr<IMnozenie> createStrassen() {
    return std::make_unique<StrassenImpl>();
}

