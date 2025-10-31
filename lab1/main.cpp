#include <iostream>
#include <chrono>
#include <iomanip>
#include <vector>
#include <memory>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <random>

#include "SupportFunctions.h"
#include "Mnozenie.h"
#include "Binet.h"   // fabryka createBinet()
#include "AI.h"      // fabryka createAI()

// create random rectangular matrix (rows x cols)
static Matrix createRandomMatrixRect(int rows, int cols, double minv = -1.0, double maxv = 1.0) {
    std::mt19937_64 rng(1234567); // deterministic seed for reproducibility
    std::uniform_real_distribution<double> dist(minv, maxv);
    Matrix M(rows, std::vector<double>(cols));
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            M[i][j] = dist(rng);
    return M;
}

// naive multiply for general rectangular matrices
static Matrix naiveMultiply(const Matrix& A, const Matrix& B) {
    int m = static_cast<int>(A.size());
    if (m == 0) return Matrix{};
    int k = static_cast<int>(A[0].size());
    int kb = static_cast<int>(B.size());
    if (k != kb) throw std::runtime_error("naiveMultiply: incompatible dims");
    int n = static_cast<int>(B[0].size());
    Matrix C = zeroMatrix(m, n);
    for (int i = 0; i < m; ++i) {
        for (int t = 0; t < k; ++t) {
            double a = A[i][t];
            for (int j = 0; j < n; ++j) {
                C[i][j] += a * B[t][j];
            }
        }
    }
    return C;
}

static std::pair<bool,double> compareMatrices(const Matrix& X, const Matrix& Y, double tol = 1e-9) {
    if (X.size() != Y.size()) return {false, std::numeric_limits<double>::infinity()};
    int m = static_cast<int>(X.size());
    if (m == 0) return {true, 0.0};
    if (X[0].size() != Y[0].size()) return {false, std::numeric_limits<double>::infinity()};
    int n = static_cast<int>(X[0].size());
    double maxDiff = 0.0;
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            double d = std::abs(X[i][j] - Y[i][j]);
            if (d > maxDiff) maxDiff = d;
            if (maxDiff > tol) return {false, maxDiff};
        }
    }
    return {true, maxDiff};
}

int main(int argc, char** argv) {
    bool batchMode = (argc >= 2);
    std::vector<int> sizes;

    int choice = 1; // default method: Binet

    if (batchMode) {
        std::string path = argv[1];
        std::ifstream in(path);
        if (!in.is_open()) {
            std::cerr << "Nie mozna otworzyc pliku: " << path << "\n";
            return 1;
        }
        std::string line;
        while (std::getline(in, line)) {
            std::istringstream iss(line);
            std::string token;
            if (!(iss >> token)) continue;
            if (token.size() > 0 && token[0] == '#') continue;
            try {
                int n = std::stoi(token);
                if (n > 0) sizes.push_back(n);
            } catch (...) {}
        }
        if (sizes.empty()) {
            std::cerr << "Brak poprawnych rozmiarow w pliku: " << path << "\n";
            return 1;
        }
        if (argc >= 3) {
            try {
                int m = std::stoi(argv[2]);
                if (m > 0) choice = m;
                else choice = 1;
            } catch (...) { choice = 1; }
        } else {
            choice = 1;
        }
    }

    if (!batchMode) {
        std::cout << "Podaj rozmiar macierzy N (NxN): ";
        int N;
        if (!(std::cin >> N) || N <= 0) {
            std::cerr << "Niepoprawny rozmiar.\n";
            return 1;
        }
        sizes.push_back(N);

        std::cout << "Wybierz metode:\n";
        std::cout << "1) Binet (rekurencyjnie bez padowania)\n";
        std::cout << "2) AI (specjalna implementacja: 4x5 * 5x5 -> 4x5)\n";
        std::cout << "Wybor (domyslnie 1): ";
        if (!(std::cin >> choice)) choice = 1;
    }

    const bool verbose = !batchMode;

    // ensure console prints use fixed notation (no scientific)
    std::cout << std::fixed << std::setprecision(6);

    if (batchMode) {
        std::ofstream out("wynik.txt");
        if (!out.is_open()) {
            std::cerr << "Nie mozna utworzyc pliku wynik.txt\n";
            return 1;
        }
        out << "# N metoda czas_s adds subs muls divs peak_bytes peak_calls correct max_diff\n";
        for (int N : sizes) {
            std::unique_ptr<IMnozenie> impl;
            std::string methodName = "Binet";
            bool useAI = (choice == 2 && N == 4);

            if (useAI) {
                impl = createAI();
                methodName = "AI";
            } else {
                impl = createBinet();
                methodName = (choice == 2) ? "Binet(fallback)" : "Binet";
            }

            Matrix A, B;
            if (useAI) {
                A = createRandomMatrixRect(4,5);
                B = createRandomMatrixRect(5,5);
            } else {
                A = createRandomMatrix(N);
                B = createRandomMatrix(N);
            }

            // dimension sanity
            if (A.empty() || B.empty() || (int)A[0].size() != (int)B.size()) {
                out << N << " " << methodName << " " << -1 << " 0 0 0 0 0 0 0 0\n";
                continue;
            }

            opCounterReset();
            memCounterReset();
            auto t0 = std::chrono::high_resolution_clock::now();
            Matrix C = impl->multiply(A, B);
            auto t1 = std::chrono::high_resolution_clock::now();
            OpCounts ops = opCounterGet();
            MemStats ms = memCounterGet();
            std::chrono::duration<double> elapsed = t1 - t0;

            Matrix Cref = naiveMultiply(A, B);
            auto [ok, maxDiff] = compareMatrices(C, Cref, 1e-9);

            // write maxDiff in fixed (no scientific)
            out << N << " " << methodName << " " << std::fixed << std::setprecision(6) << elapsed.count()
                << " " << ops.adds << " " << ops.subs << " " << ops.muls << " " << ops.divs
                << " " << ms.peak_bytes << " " << ms.peak_calls
                << " " << (ok ? "1" : "0") << " " << std::fixed << std::setprecision(6) << maxDiff << "\n";
        }
        out.close();
        return 0;
    }

    int N = sizes[0];
    bool useAI = (choice == 2 && N == 4);

    std::unique_ptr<IMnozenie> impl;
    std::string methodName = "Binet";
    if (useAI) {
        impl = createAI();
        methodName = "AI";
    } else {
        if (choice == 2 && N != 4) {
            std::cerr << "AI supports only 4x5 * 5x5. Using Binet instead.\n";
        }
        impl = createBinet();
        methodName = "Binet";
    }

    Matrix A, B;
    if (useAI) {
        A = createRandomMatrixRect(4,5);
        B = createRandomMatrixRect(5,5);
    } else {
        A = createRandomMatrix(N);
        B = createRandomMatrix(N);
    }

    opCounterReset();
    memCounterReset();
    auto t0 = std::chrono::high_resolution_clock::now();
    Matrix C = impl->multiply(A, B);
    auto t1 = std::chrono::high_resolution_clock::now();
    OpCounts ops = opCounterGet();
    MemStats ms = memCounterGet();
    std::chrono::duration<double> elapsed = t1 - t0;

    // print timing / counters to use ops/ms and avoid unused-variable warnings
    if (verbose) {
        std::cout << "Czas (s): " << std::fixed << std::setprecision(6) << elapsed.count() << "\n";
        std::cout << "Op counts: adds=" << ops.adds << " subs=" << ops.subs
                  << " muls=" << ops.muls << " divs=" << ops.divs << "\n";
        std::cout << "Memory (bytes): peak=" << ms.peak_bytes << " (peak calls=" << ms.peak_calls << ")\n";
    }

    // --- verification: compute reference result but don't print it yet ---
    Matrix Cref = naiveMultiply(A, B);
    auto [ok, maxDiff] = compareMatrices(C, Cref, 1e-9);
    // print verification summary (maxDiff) but delay printing the reference matrix until after C
    if (verbose) std::cout << "Verification: " << (ok ? "OK" : "FAIL") << ", maxDiff=" << std::fixed << std::setprecision(6) << maxDiff << "\n";

    // print matrices (works for rectangular as well)
    if (useAI) {
        if (verbose) { std::cout << "A (4x5):\n"; printSmall(A); }
        if (verbose) { std::cout << "B (5x5):\n"; printSmall(B); }
        if (verbose) { std::cout << "C = A * B (4x5):\n"; printSmall(C); }
        // print verification matrix after C
        if (verbose) { std::cout << "C_ref (verification result):\n"; printSmall(Cref); }
    } else {
        if (N <= 8) {
            if (verbose) { std::cout << "A:\n"; printSmall(A); }
            if (verbose) { std::cout << "B:\n"; printSmall(B); }
            if (verbose) { std::cout << "C = A * B:\n"; printSmall(C); }
            // print verification matrix after C
            if (verbose) { std::cout << "C_ref (verification result):\n"; printSmall(Cref); }
        } else {
            if (verbose) std::cout << "C[0][0] = " << std::fixed << std::setprecision(12) << C[0][0] << "\n";
            if (verbose) std::cout << "C[N-1][N-1] = " << std::fixed << std::setprecision(12) << C[N-1][N-1] << "\n";
            // for large matrices also print a couple of reference entries
            if (verbose) std::cout << "C_ref[0][0] = " << std::fixed << std::setprecision(12) << Cref[0][0] << "\n";
            if (verbose) std::cout << "C_ref[N-1][N-1] = " << std::fixed << std::setprecision(12) << Cref.back().back() << "\n";
        }
    }

    return 0;
}