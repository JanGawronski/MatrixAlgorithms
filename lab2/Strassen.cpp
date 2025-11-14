#include "Strassen.h"
#include "SupportFunctions.h"

#include <algorithm>
#include <stdexcept>
#include <memory>
#include <vector>

namespace {

using Matrix = ::Matrix;

Matrix multiplyRec(const Matrix &A, const Matrix &B) {
    if (rows(A) != cols(A) || cols(A) != rows(B) || rows(B) != cols(B)) {
        throw std::runtime_error("Implemented only for square matrices");
    }

    if (rows(A) == 1) {
        return A * B;
    }

    if (rows(A) % 2 == 0) {
        int size = rows(A);
        int halfSize = size / 2;
        memCounterEnterCall(size, size, 4);
    
        Matrix A11 = subMatrix(A, 0, 0, halfSize, halfSize);
        Matrix A12 = subMatrix(A, 0, halfSize, halfSize, halfSize);
        Matrix A21 = subMatrix(A, halfSize, 0, halfSize, halfSize);
        Matrix A22 = subMatrix(A, halfSize, halfSize, halfSize, halfSize);
        
        Matrix B11 = subMatrix(B, 0, 0, halfSize, halfSize);
        Matrix B12 = subMatrix(B, 0, halfSize, halfSize, halfSize);
        Matrix B21 = subMatrix(B, halfSize, 0, halfSize, halfSize);
        Matrix B22 = subMatrix(B, halfSize, halfSize, halfSize, halfSize);
        
        Matrix P1 = multiplyRec(A11 + A22, B11 + B22);
        Matrix P2 = multiplyRec(A21 + A22, B11);
        Matrix P3 = multiplyRec(A11, B12 - B22);
        Matrix P4 = multiplyRec(A22, B21 - B11);
        Matrix P5 = multiplyRec(A11 + A12, B22);
        Matrix P6 = multiplyRec(A21 - A11, B11 + B12);
        Matrix P7 = multiplyRec(A12 - A22, B21 + B22);

        Matrix C11 = P1 + P4 - P5 + P7;
        Matrix C12 = P3 + P5;
        Matrix C21 = P2 + P4;
        Matrix C22 = P1 + P3 - P2 + P6;

        Matrix M = combine(C11, C12, C21, C22);

        memCounterExitCall(size, size, 4);
        return M;
    } else {
        int size = rows(A);

        memCounterEnterCall(size, size, 4);

        Matrix A11 = subMatrix(A, 0, 0, size - 1, size - 1);
        Matrix A12 = subMatrix(A, 0, size - 1, size - 1, 1);
        Matrix A21 = subMatrix(A, size - 1, 0, 1, size - 1);
        Matrix A22 = subMatrix(A, size - 1, size - 1, 1, 1);

        Matrix B11 = subMatrix(B, 0, 0, size - 1, size - 1);
        Matrix B12 = subMatrix(B, 0, size - 1, size - 1, 1);
        Matrix B21 = subMatrix(B, size - 1, 0, 1, size - 1);
        Matrix B22 = subMatrix(B, size - 1, size - 1, 1, 1);

        Matrix C11 = multiplyRec(A11, B11) + A12 * B21;
        Matrix C12 = A11 * B12 + A12 * B22;
        Matrix C21 = A21 * B11 + A22 * B21;
        Matrix C22 = A21 * B12 + A22 * B22;

        Matrix M = combine(C11, C12, C21, C22);

        memCounterExitCall(size, size, 4);

        return M;
    }
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

