#include "GaussElimination.h"
#include "SupportFunctions.h"
#include "LUfactorization.h"
#include "Inverse.h"

std::pair<Matrix, Matrix> GaussElimination(const Matrix &A, const Matrix &b, std::unique_ptr<IMnozenie> &multImpl) {
    if (rows(A) == 1) {
        return {A, b};
    }

    if (rows(A) % 2 == 0) {
        memCounterEnterCall(rows(A), cols(A), 4);

        int halfSize = rows(A) / 2;

        Matrix A11 = subMatrix(A, 0, 0, halfSize, halfSize);
        Matrix A12 = subMatrix(A, 0, halfSize, halfSize, halfSize);
        Matrix A21 = subMatrix(A, halfSize, 0, halfSize, halfSize);
        Matrix A22 = subMatrix(A, halfSize, halfSize, halfSize, halfSize);

        Matrix b1 = subMatrix(b, 0, 0, halfSize, 1);
        Matrix b2 = subMatrix(b, halfSize, 0, halfSize, 1);

        auto [L11, U11] = LUfactorization(A11, multImpl);

        Matrix L11_inv = inverse(L11, multImpl);
        Matrix U11_inv = inverse(U11, multImpl);

        Matrix S1 = multImpl->multiply(A21, U11_inv);
        Matrix S2 = multImpl->multiply(L11_inv, A12);
        Matrix S3 = L11_inv * b1;

        auto [LS, US] = LUfactorization(A22 - multImpl->multiply(S1, S2), multImpl);

        Matrix LS_inv = inverse(LS, multImpl);

        Matrix c1 = S3;
        Matrix c2 = LS_inv * b2 - multImpl->multiply(LS_inv, S1) * S3;

        Matrix C11 = U11;
        Matrix C12 = multImpl->multiply(L11, A12);
        Matrix C21 = zeroMatrix(halfSize, halfSize);
        Matrix C22 = US;

        Matrix C = combine(C11, C12, 
                           C21, C22);
        Matrix c = combine(c1, {}, 
                           c2, {});

        memCounterExitCall(rows(A), cols(A), 4);
        return {C, c};
    } else {
        Matrix A_padded = pad(A, rows(A) + 1, rows(A) + 1);
        Matrix b_padded = pad(b, rows(b) + 1, 1);
        A_padded[rows(A)][rows(A)] = 1.0;
        b_padded[rows(b)][0] = 0.0;
        auto [C_padded, c_padded] = GaussElimination(A_padded, b_padded, multImpl);
        Matrix C = trim(C_padded, rows(A), rows(A));
        Matrix c = trim(c_padded, rows(b), 1);
        return {C, c};
    }
}