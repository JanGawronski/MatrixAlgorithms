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
    bool batchMode = (argc >= 2);
    std::vector<int> sizes;

    int choice = 1; // domyślna metoda

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
                else { choice = 1; }
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
        std::cout << "2) Strassen\n";
        std::cout << "3) AI\n";
        std::cout << "Wybor (domyslnie 1): ";
        if (!(std::cin >> choice)) choice = 1;
    }

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

    const bool verbose = !batchMode;

    if (batchMode) {
        std::ofstream out("wynik.txt");
        if (!out.is_open()) {
            std::cerr << "Nie mozna utworzyc pliku wynik.txt\n";
            return 1;
        }
        out << "# N czas_s adds subs muls divs peak_bytes peak_calls\n";
        for (int N : sizes) {
            auto A = createRandomMatrix(N);
            auto B = createRandomMatrix(N);
            if ((int)A.size() != N || (int)B.size() != N) {
                out << N << " " << -1 << " 0 0 0 0 0 0\n";
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

            out << N << " " << std::fixed << std::setprecision(6) << elapsed.count()
                << " " << ops.adds << " " << ops.subs << " " << ops.muls << " " << ops.divs
                << " " << ms.peak_bytes << " " << ms.peak_calls << "\n";
        }
        out.close();
        return 0;
    }

    int N = sizes[0];
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

    if (verbose) {
        std::cout << "Czas (s): " << std::fixed << std::setprecision(6) << elapsed.count() << "\n";
        std::cout << "Op counts: adds=" << ops.adds << " subs=" << ops.subs
                  << " muls=" << ops.muls << " divs=" << ops.divs << "\n";
        std::cout << "Memory (bytes): peak=" << ms.peak_bytes << " (peak calls=" << ms.peak_calls << ")\n";
    }

    if (verbose) {
        if (N <= 8) {
            std::cout << "A:\n"; printSmall(A);
            std::cout << "B:\n"; printSmall(B);
            std::cout << "C = A * B:\n"; printSmall(C);
            Matrix C_ref = A * B;
            Matrix diff = C - C_ref;
            std::cout << "C - C_ref:\n"; printSmall(diff);
        } else {
            std::cout << "C[0][0] = " << std::setprecision(12) << C[0][0] << "\n";
            std::cout << "C[N-1][N-1] = " << C[N-1][N-1] << "\n";
        }
    }

    return 0;
}