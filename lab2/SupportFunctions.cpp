#include "SupportFunctions.h"
#include <random>
#include <iomanip>
#include <vector>
#include <atomic>

// we use plain uint64_t (single-threaded). If multithreaded, change to atomic.
static std::uint64_t g_adds = 0;
static std::uint64_t g_subs = 0;
static std::uint64_t g_muls = 0;
static std::uint64_t g_divs = 0;

void opCounterReset() {
    g_adds = g_subs = g_muls = g_divs = 0;
}

OpCounts opCounterGet() {
    return OpCounts{g_adds, g_subs, g_muls, g_divs};
}

void opCounterAdd(const OpCounts &c) {
    g_adds += c.adds;
    g_subs += c.subs;
    g_muls += c.muls;
    g_divs += c.divs;
}

static std::uint64_t g_mem_current = 0;
static std::uint64_t g_mem_peak = 0;
static std::uint64_t g_active_calls = 0;
static std::uint64_t g_peak_calls = 0;

void memCounterReset() {
    g_mem_current = 0;
    g_mem_peak = 0;
    g_active_calls = 0;
    g_peak_calls = 0;
}

void memCounterEnterCall(int p, int r, int n) {
    std::uint64_t bytes = p * r * sizeof(double) * n;
    g_mem_current += bytes;
    if (g_mem_current > g_mem_peak) g_mem_peak = g_mem_current;
    ++g_active_calls;
    if (g_active_calls > g_peak_calls) g_peak_calls = g_active_calls;
}

void memCounterExitCall(int p, int r, int n) {
    std::uint64_t bytes = p * r * sizeof(double) * n;
    // avoid underflow
    if (g_mem_current >= bytes) g_mem_current -= bytes;
    else g_mem_current = 0;
    if (g_active_calls > 0) --g_active_calls;
}

MemStats memCounterGet() {
    return MemStats{g_mem_current, g_mem_peak, g_active_calls, g_peak_calls};
}

Matrix createRandomMatrix(int m, int n) {
    Matrix M(m, std::vector<double>(n));
    std::mt19937_64 rng(std::random_device{}());
    // open interval (1e-8, 1.0) approximated by avoiding exact endpoints
    std::uniform_real_distribution<double> dist(1e-8 + 1e-16, 1.0 - 1e-16);

    for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j)
            M[i][j] = dist(rng);
    return M;
}

Matrix createRandomMatrix(int n) {
    return createRandomMatrix(n, n);
}


void printSmall(const Matrix& M) {
    for (size_t i = 0; i < M.size(); ++i) {
        for (size_t j = 0; j < M[i].size(); ++j)
            std::cout << std::setprecision(6) << std::setw(12) << M[i][j];
        std::cout << '\n';
    }
}

Matrix zeroMatrix(int rows, int cols) {
    return Matrix(rows, std::vector<double>(cols, 0.0));
}

Matrix zeroMatrix(int n) {
    return zeroMatrix(n, n);
}

// A[row:row+rows-1][col:col+cols-1]
Matrix subMatrix(const Matrix &A, int row, int col, int rows, int cols) {
    Matrix R = zeroMatrix(rows, cols);
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            R[i][j] = A[row + i][col + j];
    return R;
}

Matrix combine(const Matrix &A11, const Matrix &A12,
               const Matrix &A21, const Matrix &A22) {
    Matrix A = zeroMatrix(rows(A11) + rows(A21), cols(A11) + cols(A12));
    for (int i = 0; i < rows(A11); i++)
        for (int j = 0; j < cols(A11); j++)
            A[i][j] = A11[i][j];

    for (int i = 0; i < rows(A12); i++)
        for (int j = 0; j < cols(A12); j++)
            A[i][j + cols(A11)] = A12[i][j];

    for (int i = 0; i < rows(A21); i++)
        for (int j = 0; j < cols(A11); j++)
            A[i + rows(A11)][j] = A21[i][j];

    for (int i = 0; i < rows(A22); i++)
        for (int j = 0; j < cols(A12); j++)
            A[i + rows(A11)][j + cols(A11)] = A22[i][j];

    return A;
}

int rows(const Matrix& M) {
    return static_cast<int>(M.size());
}

int cols(const Matrix& M) {
    return M.size() ? static_cast<int>(M[0].size()) : 0;
}

Matrix operator+(const Matrix &A, const Matrix &B) {
    Matrix R = zeroMatrix(std::max(rows(A), rows(B)), std::max(cols(A), cols(B)));
    for (int i = 0; i < rows(A); ++i)
        for (int j = 0; j < cols(A); ++j)
            R[i][j] = A[i][j];

    for (int i = 0; i < rows(B); ++i) {
        for (int j = 0; j < cols(B); ++j) {
            R[i][j] += B[i][j];
            ++g_adds;
        }
    }

    return R;
}

Matrix operator-(const Matrix &A, const Matrix &B) {
    Matrix R = zeroMatrix(std::max(rows(A), rows(B)), std::max(cols(A), cols(B)));
    for (int i = 0; i < rows(A); ++i)
        for (int j = 0; j < cols(A); ++j)
            R[i][j] = A[i][j];

    for (int i = 0; i < rows(B); ++i) {
        for (int j = 0; j < cols(B); ++j) {
            R[i][j] -= B[i][j];
            ++g_subs;
        }
    }

    return R;
}

Matrix operator*(const Matrix &A, const Matrix &B) {
    memCounterEnterCall(rows(A), cols(B), 1);

    Matrix C = zeroMatrix(rows(A), cols(B));
    for (int i = 0; i < rows(A); ++i) {
        for (int k = 0; k < cols(A); ++k) {
            double aik = A[i][k];
            for (int j = 0; j < cols(B); ++j) {
                double prod = aik * B[k][j];
                ++g_muls;
                C[i][j] += prod;
                ++g_adds;
            }
        }
    }
    g_adds -= rows(A) * cols(B);

    memCounterExitCall(rows(A), cols(B), 1);
    return C;
}

std::pair<bool,double> compareMatrices(const Matrix& X, const Matrix& Y, double tol) {
    if (rows(X) != rows(Y)) return {false, std::numeric_limits<double>::infinity()};
    if (rows(X) == 0) return {true, 0.0};
    if (X[0].size() != Y[0].size()) return {false, std::numeric_limits<double>::infinity()};
    int n = static_cast<int>(X[0].size());
    double maxDiff = 0.0;
    for (int i = 0; i < rows(X); ++i) {
        for (int j = 0; j < n; ++j) {
            double d = std::abs(X[i][j] - Y[i][j]);
            if (d > maxDiff) maxDiff = d;
            if (maxDiff > tol) return {false, maxDiff};
        }
    }
    return {true, maxDiff};
}
