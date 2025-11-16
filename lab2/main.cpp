#include <iostream>
#include <chrono>
#include <iomanip>
#include <vector>
#include <memory>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

#include "SupportFunctions.h"
#include "Mnozenie.h"
#include "Binet.h"
#include "Strassen.h"
#include "Inverse.h"
#include "LUfactorization.h"
#include "GaussElimination.h"

int main(int argc, char** argv) {
    if (argc >= 2) { //batch mode
        std::vector<int> sizes;
        std::string path;
        if (argc >= 3) {
            path = argv[2];
        } else {
            path = "sizes.txt";
        }
        std::ifstream in(path);
        if (!in.is_open()) {
            std::cerr << "Can't open file: " << path << "\n";
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
        in.close();
        
        if (sizes.empty()) {
            std::cerr << "No correct sizes: " << path << "\n";
            return 1;
        }

        std::string outpath;
        if (argc >= 4) {
            outpath = argv[3];
        } else {
            outpath = "results.txt";
        }

        int choice = argv[1][0] - '0';

        std::unique_ptr<IMnozenie> impl = choice % 2 == 1 ? createBinet() : createStrassen();

        std::ofstream out(outpath);
        if (!out.is_open()) {
            std::cerr << "Can't open results file\n";
            return 1;
        }

        out << "# N time_s adds subs muls divs peak_bytes peak_calls\n";

        for (size_t i = 0; i < sizes.size(); ++i) {
            int N = sizes[i];
            
            std::cout << "[" << (i+1) << "/" << sizes.size() << "] ";
            std::cout << "N=" << std::setw(4) << N << " ... ";
            std::cout.flush();

            Matrix A = createRandomMatrix(N);
            Matrix b = createRandomMatrix(N, 1);

            opCounterReset();
            memCounterReset();
            auto t0 = std::chrono::high_resolution_clock::now();

            switch ((choice - 1) / 2)
            {
            case 0:
                inverse(A, impl);
                break;
            case 1:
                LUfactorization(A, impl);
                break;
            case 2:
                GaussElimination(A, b, impl);
                break;
            default:
                break;
            }
            
            auto t1 = std::chrono::high_resolution_clock::now();
            OpCounts ops = opCounterGet();
            MemStats ms = memCounterGet();
            std::chrono::duration<double> elapsed = t1 - t0;

            out << N << " " << std::fixed << std::setprecision(6) << elapsed.count()
                << " " << ops.adds << " " << ops.subs << " " << ops.muls << " " << ops.divs
                << " " << ms.peak_bytes << " " << ms.peak_calls << "\n";
            out.flush();

            std::cout << std::fixed << std::setprecision(3) << elapsed.count() << "s | ";
            std::cout.flush();
            
            std::cout << " OK\n";
            std::cout.flush();
        }

        out.close();
    } else {
        std::cout << "Enter matrix size N (NxN): ";
        int N;
        if (!(std::cin >> N) || N <= 0) {
            std::cerr << "Incorrect size.\n";
            return 1;
        }

        int choice;
        std::cout << "Choose operation:\n";
        std::cout << "1) Inverse matrix (Binet)\n";
        std::cout << "2) Inverse matrix (Strassen)\n";
        std::cout << "3) Gauss elimination (Binet)\n";
        std::cout << "4) Gauss elimination (Strassen)\n";
        std::cout << "5) LU factorization (Binet)\n";
        std::cout << "6) LU factorization (Strassen)\n";
        std::cout << "Choice: ";
        if (!(std::cin >> choice)) choice = 1;

        Matrix A = createRandomMatrix(N);
        std::unique_ptr<IMnozenie> impl = choice % 2 == 1 ? createBinet() : createStrassen();

        switch ((choice - 1) / 2) {
            case 0: {
                opCounterReset();
                memCounterReset();
                auto t0 = std::chrono::high_resolution_clock::now();
                Matrix B = inverse(A, impl);
                auto t1 = std::chrono::high_resolution_clock::now();
                OpCounts ops = opCounterGet();
                MemStats ms = memCounterGet();
                std::chrono::duration<double> elapsed = t1 - t0;

                std::cout << "Time (s): " << std::fixed << std::setprecision(6) << elapsed.count() << "\n";
                std::cout << "Op counts: adds=" << ops.adds << " subs=" << ops.subs
                            << " muls=" << ops.muls << " divs=" << ops.divs << "\n";
                std::cout << "Memory (bytes): peak=" << ms.peak_bytes << " (peak calls=" << ms.peak_calls << ")\n";
                
                B = inverse(A, impl);
                Matrix C = A * B;
                auto [equal, max_err] = compareMatrices(C, identityMatrix(N), 1e-6);
                if (equal) {
                    std::cout << "Inverse check passed (max error=" << std::setprecision(3) << max_err << ")\n";
                } else {
                    std::cout << "Inverse check FAILED (max error=" << std::setprecision(3) << max_err << ")\n";
                }
                if (N <= 12) {
                    std::cout << "A:\n"; printSmall(A);
                    std::cout << "B:\n"; printSmall(B);
                    std::cout << "A*B:\n"; printSmall(C);
                } else {
                    std::cout << "B[0][0] = " << std::setprecision(12) << B[0][0] << "\n";
                    std::cout << "B[N-1][N-1] = " << B[N-1][N-1] << "\n";
                }
                break;
            }
            case 1: {
                Matrix b = createRandomMatrix(N, 1);
                opCounterReset();
                memCounterReset();
                auto t0 = std::chrono::high_resolution_clock::now();
                auto [C, c] = GaussElimination(A, b, impl);
                auto t1 = std::chrono::high_resolution_clock::now();
                OpCounts ops = opCounterGet();
                MemStats ms = memCounterGet();
                std::chrono::duration<double> elapsed = t1 - t0;

                std::cout << "Time (s): " << std::fixed << std::setprecision(6) << elapsed.count() << "\n";
                std::cout << "Op counts: adds=" << ops.adds << " subs=" << ops.subs
                            << " muls=" << ops.muls << " divs=" << ops.divs << "\n";
                std::cout << "Memory (bytes): peak=" << ms.peak_bytes << " (peak calls=" << ms.peak_calls << ")\n";

                if (N <= 12) {
                    std::cout << "A:\n"; printSmall(A);
                    std::cout << "b:\n"; printSmall(b);
                    std::cout << "C:\n"; printSmall(C);
                    std::cout << "c:\n"; printSmall(c);
                } else {
                    std::cout << "A[0][0] = " << std::setprecision(12) << A[0][0] << "\n";
                    std::cout << "A[N-1][N-1] = " << A[N-1][N-1] << "\n";
                }
                break;
            }
            case 2: {
                opCounterReset();
                memCounterReset();
                auto t0 = std::chrono::high_resolution_clock::now();
                auto [L, U] = LUfactorization(A, impl);
                auto t1 = std::chrono::high_resolution_clock::now();
                OpCounts ops = opCounterGet();
                MemStats ms = memCounterGet();
                std::chrono::duration<double> elapsed = t1 - t0;

                std::cout << "Time (s): " << std::fixed << std::setprecision(6) << elapsed.count() << "\n";
                std::cout << "Op counts: adds=" << ops.adds << " subs=" << ops.subs
                            << " muls=" << ops.muls << " divs=" << ops.divs << "\n";
                std::cout << "Memory (bytes): peak=" << ms.peak_bytes << " (peak calls=" << ms.peak_calls << ")\n";
                std::cout << "Determinant: " << determinantLU(A, impl) << "\n";

                Matrix LU = impl->multiply(L, U);
                auto [equal, max_err] = compareMatrices(LU, A, 1e-6);
                if (equal) {
                    std::cout << "LU factorization check passed (max error=" << std::setprecision(3) << max_err << ")\n";
                } else {
                    std::cout << "LU factorization check FAILED (max error=" << std::setprecision(3) << max_err << ")\n";
                }
                if (N <= 12) {
                    std::cout << "A:\n"; printSmall(A);
                    std::cout << "L:\n"; printSmall(L);
                    std::cout << "U:\n"; printSmall(U);
                    std::cout << "L*U:\n"; printSmall(LU);
                    std::cout << "L*U - A:\n"; printSmall(LU - A);
                } else {
                    std::cout << "L[0][0] = " << std::setprecision(12) << L[0][0] << "\n";
                    std::cout << "U[N-1][N-1] = " << U[N-1][N-1] << "\n";
                }
                break;
            }
            default:
                std::cerr << "Incorrect method.\n";
                return 1;
        }

    }

    return 0;
}