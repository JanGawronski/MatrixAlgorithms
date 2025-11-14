#include "Binet.h"
#include "SupportFunctions.h"

#include <algorithm>
#include <stdexcept>
#include <memory>
#include <vector>

namespace {

using Matrix = ::Matrix;

Matrix multiplyRec(const Matrix &A, const Matrix &B) {
    if (cols(A) != rows(B)) {
        throw std::runtime_error("Incompatible dimensions for multiplication");
    }

    if (rows(A) == 1 || cols(A) == 1 || cols(B) == 1) {
        Matrix M = A * B;
        return M;
    }

    memCounterEnterCall(rows(A), cols(B), 3);

    int A11width = cols(A) / 2;
    int A12width = cols(A) - A11width;
    int A11height = rows(A) / 2;
    int A21height = rows(A) - A11height;

    int B11width = cols(B) / 2;
    int B12width = cols(B) - B11width;
    
    Matrix A11 = subMatrix(A, 0, 0, A11height, A11width);
    Matrix A12 = subMatrix(A, 0, A11width, A11height, A12width);
    Matrix A21 = subMatrix(A, A11height, 0, A21height, A11width);
    Matrix A22 = subMatrix(A, A11height, A11width, A21height, A12width);

    Matrix B11 = subMatrix(B, 0, 0, A11width, B11width);
    Matrix B12 = subMatrix(B, 0, B11width, A11width, B12width);
    Matrix B21 = subMatrix(B, A11width, 0, A12width, B11width);
    Matrix B22 = subMatrix(B, A11width, B11width, A12width, B12width);

    Matrix C11 = multiplyRec(A11, B11) + multiplyRec(A12, B21);
    Matrix C12 = multiplyRec(A11, B12) + multiplyRec(A12, B22);
    Matrix C21 = multiplyRec(A21, B11) + multiplyRec(A22, B21);
    Matrix C22 = multiplyRec(A21, B12) + multiplyRec(A22, B22);

    Matrix M = combine(C11, C12, C21, C22);

    memCounterExitCall(rows(A), cols(B), 3);
    return M;
}

} // namespace (internal)


class BinetImpl : public IMnozenie {
public:
    Matrix multiply(const Matrix &A, const Matrix &B) override {
        return multiplyRec(A, B);
    }
};

std::unique_ptr<IMnozenie> createBinet() {
    return std::make_unique<BinetImpl>();
}

