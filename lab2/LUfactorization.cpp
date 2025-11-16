#include "LUfactorization.h"
#include "SupportFunctions.h"
#include "Inverse.h"

std::pair<Matrix, Matrix> LUfactorization(const Matrix &A, std::unique_ptr<IMnozenie> &multImpl) {
    if (rows(A) == 1) {
        Matrix L = identityMatrix(1);
        Matrix U = A;
        return {L, U};
    }
    
    if (rows(A) % 2 == 0) {
        memCounterEnterCall(rows(A), cols(A), 3);

        int halfSize = rows(A) / 2;

        Matrix A11 = subMatrix(A, 0, 0, halfSize, halfSize);
        Matrix A12 = subMatrix(A, 0, halfSize, halfSize, halfSize);
        Matrix A21 = subMatrix(A, halfSize, 0, halfSize, halfSize);
        Matrix A22 = subMatrix(A, halfSize, halfSize, halfSize, halfSize);

        auto [L11, U11] = LUfactorization(A11, multImpl);

        Matrix L11_inv = inverse(L11, multImpl);
        Matrix U11_inv = inverse(U11, multImpl);

        Matrix U12 = multImpl->multiply(L11_inv, A12);
        Matrix L21 = multImpl->multiply(A21, U11_inv);

        Matrix S = A22 - multImpl->multiply(L21, U12);

        auto [L22, U22] = LUfactorization(S, multImpl);

        Matrix L = combine(L11, zeroMatrix(halfSize, halfSize),
                           L21, L22);
        Matrix U = combine(U11, U12,
                           zeroMatrix(halfSize, halfSize), U22);

        memCounterExitCall(rows(A), cols(A), 3);
        return {L, U};
    } else {
        Matrix A_padded = pad(A, rows(A) + 1, rows(A) + 1);
        auto [L_padded, U_padded] = LUfactorization(A_padded, multImpl);
        Matrix L = trim(L_padded, rows(A), rows(A));
        Matrix U = trim(U_padded, rows(A), rows(A));
        return {L, U};
    }
}

double determinantLU(const Matrix &A, std::unique_ptr<IMnozenie> &multImpl) {
    auto [_, U] = LUfactorization(A, multImpl);
    double det = 1.0;
    for (int i = 0; i < rows(U); ++i) {
        det *= U[i][i];
    }
    return det;
}