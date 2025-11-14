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

        std::unique_ptr<IMnozenie> multiplyImpl;
        switch (argv[1][0]) {
            case '1':
                multiplyImpl = createBinet();
                break;
            case '2':
                multiplyImpl = createStrassen();
                break;
            default:
                std::cerr << "Unknown method: " << argv[1] << "\n";
                return 1;
        }

        std::ofstream out(outpath);
        if (!out.is_open()) {
            std::cerr << "Can't open results file\n";
            return 1;
        }

        out << "# N czas_s adds subs muls divs peak_bytes peak_calls\n";

        for (size_t i = 0; i < sizes.size(); ++i) {
            int N = sizes[i];
            
            std::cout << "[" << (i+1) << "/" << sizes.size() << "] ";
            std::cout << "N=" << std::setw(4) << N << " ... ";
            std::cout.flush();

            auto A = createRandomMatrix(N);
            auto B = createRandomMatrix(N);
            if (rows(A) != N || rows(B) != N) {
                out << N << " " << -1 << " 0 0 0 0 0 0\n";
                std::cout << "ERROR (creating matrix)\n";
                continue;
            }

            
            opCounterReset();
            memCounterReset();
            auto t0 = std::chrono::high_resolution_clock::now();
            Matrix C = multiplyImpl->multiply(A, B);
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