#include "Binet.h"
#include "SupportFunctions.h"

#include <algorithm>
#include <stdexcept>
#include <memory>
#include <vector>

// wszystkie pomocnicze funkcje i symbole o zasięgu pliku
namespace {

using Matrix = ::Matrix; // użyj aliasu zdefiniowanego w Mnozenie.h

Matrix combineVertical(const Matrix &top, const Matrix &bottom) {
    int topR = static_cast<int>(top.size()), botR = static_cast<int>(bottom.size());
    int cols = topR ? static_cast<int>(top[0].size()) : (botR ? static_cast<int>(bottom[0].size()) : 0);
    Matrix R = zeroMatrix(topR + botR, cols);
    for (int i = 0; i < topR; ++i)
        for (int j = 0; j < cols; ++j)
            R[i][j] = top[i][j];
    for (int i = 0; i < botR; ++i)
        for (int j = 0; j < cols; ++j)
            R[topR + i][j] = bottom[i][j];
    return R;
}

Matrix combineHorizontal(const Matrix &left, const Matrix &right) {
    int rows = static_cast<int>(left.size());
    int leftC = rows ? static_cast<int>(left[0].size()) : 0;
    int rightC = rows ? static_cast<int>(right[0].size()) : 0;
    Matrix R = zeroMatrix(rows, leftC + rightC);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < leftC; ++j) R[i][j] = left[i][j];
        for (int j = 0; j < rightC; ++j) R[i][leftC + j] = right[i][j];
    }
    return R;
}

// rekurencyjne mnożenie (bez pad'owania) - dzieli po największym wymiarze
Matrix multiplyRec(const Matrix &A, const Matrix &B) {
    int p = static_cast<int>(A.size());
    int q = p ? static_cast<int>(A[0].size()) : 0; // A is p x q
    int qb = static_cast<int>(B.size());           // B is qb x r
    int r = qb ? static_cast<int>(B[0].size()) : 0;
    if (q != qb) {
        throw std::runtime_error("Incompatible dimensions for multiplication");
    }

    // account memory for this call (result p x r)
    memCounterEnterCall(static_cast<std::size_t>(p), static_cast<std::size_t>(r));

    const int BASE = 4;
    if (p == 0 || q == 0 || r == 0) {
        Matrix Z = zeroMatrix(p, r);
        memCounterExitCall(static_cast<std::size_t>(p), static_cast<std::size_t>(r));
        return Z;
    }
    if (std::max({p, q, r}) <= BASE) {
        Matrix res = multiplyClassic(A, B);
        memCounterExitCall(static_cast<std::size_t>(p), static_cast<std::size_t>(r));
        return res;
    }

    if (p >= std::max(q, r)) {
        int p1 = p / 2;
        int p2 = p - p1;
        Matrix A_top = subMatrix(A, 0, 0, p1, q);
        Matrix A_bot = subMatrix(A, p1, 0, p2, q);
        Matrix C_top = multiplyRec(A_top, B);
        Matrix C_bot = multiplyRec(A_bot, B);
        Matrix R = combineVertical(C_top, C_bot);
        memCounterExitCall(static_cast<std::size_t>(p), static_cast<std::size_t>(r));
        return R;
    } else if (r >= std::max(p, q)) {
        int r1 = r / 2;
        int r2 = r - r1;
        Matrix B_left = subMatrix(B, 0, 0, q, r1);
        Matrix B_right = subMatrix(B, 0, r1, q, r2);
        Matrix C_left = multiplyRec(A, B_left);
        Matrix C_right = multiplyRec(A, B_right);
        Matrix R = combineHorizontal(C_left, C_right);
        memCounterExitCall(static_cast<std::size_t>(p), static_cast<std::size_t>(r));
        return R;
    } else {
        int q1 = q / 2;
        int q2 = q - q1;
        Matrix A_left = subMatrix(A, 0, 0, p, q1);
        Matrix A_right = subMatrix(A, 0, q1, p, q2);
        Matrix B_top = subMatrix(B, 0, 0, q1, r);
        Matrix B_bot = subMatrix(B, q1, 0, q2, r);
        Matrix P1 = multiplyRec(A_left, B_top);
        Matrix P2 = multiplyRec(A_right, B_bot);
        Matrix R = addMatrix(P1, P2);
        memCounterExitCall(static_cast<std::size_t>(p), static_cast<std::size_t>(r));
        return R;
    }
}

} // namespace (internal)


// IMnozenie implementation hidden in this translation unit
class BinetImpl : public IMnozenie {
public:
    Matrix multiply(const Matrix &A, const Matrix &B) override {
        return multiplyRec(A, B);
    }
};

std::unique_ptr<IMnozenie> createBinet() {
    return std::make_unique<BinetImpl>();
}

