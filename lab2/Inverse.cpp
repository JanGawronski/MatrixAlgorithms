#include "Inverse.h"
#include "SupportFunctions.h"

Matrix inverse(const Matrix &A, std::unique_ptr<IMnozenie> &multImpl) {
    if (rows(A) == 1) {
        Matrix invA = zeroMatrix(1, 1);
        invA[0][0] = 1.0 / A[0][0];
        opCounterAdd({0, 0, 0, 1});
        return invA;
    }

    if (rows(A) % 2 == 0) {

        memCounterEnterCall(rows(A), cols(A), 3);

        int halfSize = rows(A) / 2;
        Matrix invA11 = inverse(subMatrix(A, 0, 0, halfSize, halfSize), multImpl);
        Matrix A12 = subMatrix(A, 0, halfSize, halfSize, halfSize);
        Matrix A21 = subMatrix(A, halfSize, 0, halfSize, halfSize);
        Matrix A22 = subMatrix(A, halfSize, halfSize, halfSize, halfSize);

        Matrix T1 = multImpl->multiply(invA11, A12);
        Matrix T2 = multImpl->multiply(A21, invA11);
        
        Matrix invS22 = inverse(A22 - multImpl->multiply(A21, T1), multImpl);

        Matrix T3 = multImpl->multiply(T1, invS22);

        Matrix B11 = invA11 + multImpl->multiply(T3, T2);
        Matrix B12 = negate(T3);
        Matrix B21 = negate(multImpl->multiply(invS22, T2));
        Matrix B22 = invS22;

        memCounterExitCall(rows(A), cols(A), 3);
        return combine(B11, B12, B21, B22);
    } else {
        Matrix A_padded = pad(A, rows(A) + 1, cols(A) + 1);
        A_padded[rows(A)][cols(A)] = 1.0;
        Matrix inv_padded = inverse(A_padded, multImpl);
        return trim(inv_padded, rows(A), cols(A));
    }
}