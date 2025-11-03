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
#include "AI.h"

int main(int argc, char** argv) {
    if (argc >= 2) { //batch mode
        std::vector<int> sizes;
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

        // Create both implementations
        auto binetImpl = createBinet();
        auto strassenImpl = createStrassen();

        // Open output files
        std::ofstream outBinet("wynikBinet.txt");
        std::ofstream outStrassen("wynikStrassen.txt");
        if (!outBinet.is_open() || !outStrassen.is_open()) {
            std::cerr << "Nie mozna utworzyc plikow wynikowych\n";
            return 1;
        }

        outBinet << "# N czas_s adds subs muls divs peak_bytes peak_calls\n";
        outStrassen << "# N czas_s adds subs muls divs peak_bytes peak_calls\n";

        std::cout << "\n========================================\n";
        std::cout << "Rozpoczynam obliczenia dla " << sizes.size() << " rozmiarow macierzy\n";
        std::cout << "========================================\n\n";

        // Process each size alternating between Binet and Strassen
        for (size_t i = 0; i < sizes.size(); ++i) {
            int N = sizes[i];
            
            std::cout << "[" << (i+1) << "/" << sizes.size() << "] ";
            std::cout << "N=" << std::setw(4) << N << " ... ";
            std::cout.flush();

            auto A = createRandomMatrix(N);
            auto B = createRandomMatrix(N);
            if ((int)A.size() != N || (int)B.size() != N) {
                outBinet << N << " " << -1 << " 0 0 0 0 0 0\n";
                outStrassen << N << " " << -1 << " 0 0 0 0 0 0\n";
                std::cout << "BLAD (tworzenie macierzy)\n";
                continue;
            }

            // Alternate: even indices -> Binet first, odd indices -> Strassen first
            bool binetFirst = (i % 2 == 0);

            if (binetFirst) {
                // Run Binet
                std::cout << "Binet... ";
                std::cout.flush();
                
                opCounterReset();
                memCounterReset();
                auto t0 = std::chrono::high_resolution_clock::now();
                Matrix C = binetImpl->multiply(A, B);
                auto t1 = std::chrono::high_resolution_clock::now();
                OpCounts ops = opCounterGet();
                MemStats ms = memCounterGet();
                std::chrono::duration<double> elapsed = t1 - t0;

                outBinet << N << " " << std::fixed << std::setprecision(6) << elapsed.count()
                    << " " << ops.adds << " " << ops.subs << " " << ops.muls << " " << ops.divs
                    << " " << ms.peak_bytes << " " << ms.peak_calls << "\n";
                outBinet.flush();

                std::cout << std::fixed << std::setprecision(3) << elapsed.count() << "s | ";
                std::cout.flush();

                // Run Strassen
                std::cout << "Strassen... ";
                std::cout.flush();
                
                opCounterReset();
                memCounterReset();
                t0 = std::chrono::high_resolution_clock::now();
                C = strassenImpl->multiply(A, B);
                t1 = std::chrono::high_resolution_clock::now();
                ops = opCounterGet();
                ms = memCounterGet();
                elapsed = t1 - t0;

                outStrassen << N << " " << std::fixed << std::setprecision(6) << elapsed.count()
                    << " " << ops.adds << " " << ops.subs << " " << ops.muls << " " << ops.divs
                    << " " << ms.peak_bytes << " " << ms.peak_calls << "\n";
                outStrassen.flush();

                std::cout << std::fixed << std::setprecision(3) << elapsed.count() << "s";
                
            } else {
                // Run Strassen first
                std::cout << "Strassen... ";
                std::cout.flush();
                
                opCounterReset();
                memCounterReset();
                auto t0 = std::chrono::high_resolution_clock::now();
                Matrix C = strassenImpl->multiply(A, B);
                auto t1 = std::chrono::high_resolution_clock::now();
                OpCounts ops = opCounterGet();
                MemStats ms = memCounterGet();
                std::chrono::duration<double> elapsed = t1 - t0;

                outStrassen << N << " " << std::fixed << std::setprecision(6) << elapsed.count()
                    << " " << ops.adds << " " << ops.subs << " " << ops.muls << " " << ops.divs
                    << " " << ms.peak_bytes << " " << ms.peak_calls << "\n";
                outStrassen.flush();

                std::cout << std::fixed << std::setprecision(3) << elapsed.count() << "s | ";
                std::cout.flush();

                // Run Binet
                std::cout << "Binet... ";
                std::cout.flush();
                
                opCounterReset();
                memCounterReset();
                t0 = std::chrono::high_resolution_clock::now();
                C = binetImpl->multiply(A, B);
                t1 = std::chrono::high_resolution_clock::now();
                ops = opCounterGet();
                ms = memCounterGet();
                elapsed = t1 - t0;

                outBinet << N << " " << std::fixed << std::setprecision(6) << elapsed.count()
                    << " " << ops.adds << " " << ops.subs << " " << ops.muls << " " << ops.divs
                    << " " << ms.peak_bytes << " " << ms.peak_calls << "\n";
                outBinet.flush();

                std::cout << std::fixed << std::setprecision(3) << elapsed.count() << "s";
            }
            
            std::cout << " OK\n";
            std::cout.flush();
        }

        outBinet.close();
        outStrassen.close();
        
        std::cout << "\n========================================\n";
        std::cout << "Zakonczono! Wyniki zapisane do:\n";
        std::cout << "  - wynikBinet.txt\n";
        std::cout << "  - wynikStrassen.txt\n";
        std::cout << "========================================\n";

    } else {
        std::cout << "Podaj rozmiar macierzy N (NxN): ";
        int N;
        if (!(std::cin >> N) || N <= 0) {
            std::cerr << "Niepoprawny rozmiar.\n";
            return 1;
        }

        int choice;
        std::cout << "Wybierz metode:\n";
        std::cout << "1) Binet (rekurencyjnie bez padowania)\n";
        std::cout << "2) Strassen\n";
        std::cout << "3) AI\n";
        std::cout << "Wybor (domyslnie 1): ";
        if (!(std::cin >> choice)) choice = 1;

        std::unique_ptr<IMnozenie> impl;
        switch (choice) {
            case 1:
                impl = createBinet();
                break;
            case 2:
                impl = createStrassen();
                break;
            case 3:
                impl = createAI();
                break;
            default:
                std::cerr << "Wybrana metoda (" << choice << ") niezaimplementowana. Uzywam Binet (1).\n";
                impl = createBinet();
        }
    
        auto A = createRandomMatrix(N);
        auto B = createRandomMatrix(N);

        if (choice == 3) {
            A = createRandomMatrix(4, 5);
            B = createRandomMatrix(5, 5);        
        }

        opCounterReset();
        memCounterReset();
        auto t0 = std::chrono::high_resolution_clock::now();
        Matrix C = impl->multiply(A, B);
        auto t1 = std::chrono::high_resolution_clock::now();
        OpCounts ops = opCounterGet();
        MemStats ms = memCounterGet();
        std::chrono::duration<double> elapsed = t1 - t0;

        std::cout << "Czas (s): " << std::fixed << std::setprecision(6) << elapsed.count() << "\n";
        std::cout << "Op counts: adds=" << ops.adds << " subs=" << ops.subs
                    << " muls=" << ops.muls << " divs=" << ops.divs << "\n";
        std::cout << "Memory (bytes): peak=" << ms.peak_bytes << " (peak calls=" << ms.peak_calls << ")\n";

        Matrix C_ref = A * B;
        auto [ok, maxdiff] = compareMatrices(C, C_ref, 1e-9);
        if (ok) {
            std::cout << "Macierze zgodne (roznica max " << std::setprecision(12) << maxdiff << ")\n";
        } else {
            std::cout << "Macierze NIEZGODNE! (roznica max " << std::setprecision(12) << maxdiff << ")\n";
        }
        
        if (N <= 12) {
            std::cout << "A:\n"; printSmall(A);
            std::cout << "B:\n"; printSmall(B);
            std::cout << "C = A * B:\n"; printSmall(C);
            std::cout << "C - C_ref:\n"; printSmall(C - C_ref);
        } else {
            std::cout << "C[0][0] = " << std::setprecision(12) << C[0][0] << "\n";
            std::cout << "C[N-1][N-1] = " << C[N-1][N-1] << "\n";
        }
    }

    return 0;
}