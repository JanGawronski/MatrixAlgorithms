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

void memCounterEnterCall(std::size_t p, std::size_t r) {
    std::uint64_t bytes = static_cast<std::uint64_t>(p) * static_cast<std::uint64_t>(r) * sizeof(double);
    g_mem_current += bytes;
    if (g_mem_current > g_mem_peak) g_mem_peak = g_mem_current;
    ++g_active_calls;
    if (g_active_calls > g_peak_calls) g_peak_calls = g_active_calls;
}

void memCounterExitCall(std::size_t p, std::size_t r) {
    std::uint64_t bytes = static_cast<std::uint64_t>(p) * static_cast<std::uint64_t>(r) * sizeof(double);
    // avoid underflow
    if (g_mem_current >= bytes) g_mem_current -= bytes;
    else g_mem_current = 0;
    if (g_active_calls > 0) --g_active_calls;
}

MemStats memCounterGet() {
    return MemStats{g_mem_current, g_mem_peak, g_active_calls, g_peak_calls};
}

Matrix createRandomMatrix(int n) {
    Matrix M(n, std::vector<double>(n));
    std::mt19937_64 rng(std::random_device{}());
    // open interval (1e-8, 1.0) approximated by avoiding exact endpoints
    std::uniform_real_distribution<double> dist(1e-8 + 1e-16, 1.0 - 1e-16);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            M[i][j] = dist(rng);
    return M;
}

void printSmall(const Matrix& M) {
    for (size_t i = 0; i < M.size(); ++i) {
        for (size_t j = 0; j < M[i].size(); ++j)
            std::cout << std::setw(12) << M[i][j];
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

Matrix operator+(const Matrix &A, const Matrix &B) {
    int rows = static_cast<int>(A.size());
    int cols = rows ? static_cast<int>(A[0].size()) : 0;
    Matrix R = zeroMatrix(rows, cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            R[i][j] = A[i][j] + B[i][j];
            ++g_adds;
        }
    }
    return R;
}

Matrix operator-(const Matrix &A, const Matrix &B) {
    int rows = static_cast<int>(A.size());
    int cols = rows ? static_cast<int>(A[0].size()) : 0;
    Matrix R = zeroMatrix(rows, cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            R[i][j] = A[i][j] - B[i][j];
            ++g_subs;
        }
    }
    return R;
}

Matrix operator*(const Matrix &A, const Matrix &B) {
    int p = static_cast<int>(A.size());
    int q = p ? static_cast<int>(A[0].size()) : 0;
    int r = B.size() ? static_cast<int>(B[0].size()) : 0;

    // account for this call's result allocation
    memCounterEnterCall(static_cast<std::size_t>(p), static_cast<std::size_t>(r));

    Matrix C = zeroMatrix(p, r);
    for (int i = 0; i < p; ++i) {
        for (int k = 0; k < q; ++k) {
            double aik = A[i][k];
            for (int j = 0; j < r; ++j) {
                double prod = aik * B[k][j];
                ++g_muls;
                C[i][j] += prod;
                ++g_adds;
            }
        }
    }

    memCounterExitCall(static_cast<std::size_t>(p), static_cast<std::size_t>(r));
    return C;
}