#include "HMatrix.h"
#include <iostream>
#include <iomanip>
#include <cmath>

int main() {
    std::cout << "=== Mnożenie H-macierzy przez wektor [1,2,3,4,5,6,7,8] ===" << std::endl;
    std::cout << std::endl;
    
    // Zbuduj przykładową macierz 8x8
    // Możesz tu wpisać konkretną macierz z tablicy
    Matrix A = {
        {4, 1, 0, 0, 1, 0, 0, 0},
        {1, 4, 1, 0, 0, 1, 0, 0},
        {0, 1, 4, 1, 0, 0, 1, 0},
        {0, 0, 1, 4, 0, 0, 0, 1},
        {1, 0, 0, 0, 4, 1, 0, 0},
        {0, 1, 0, 0, 1, 4, 1, 0},
        {0, 0, 1, 0, 0, 1, 4, 1},
        {0, 0, 0, 1, 0, 0, 1, 4}
    };
    
    // Wektor [1,2,3,4,5,6,7,8]
    Vector x = {1, 2, 3, 4, 5, 6, 7, 8};
    
    std::cout << "Macierz A (8x8):" << std::endl;
    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < 8; ++j) {
            std::cout << std::setw(4) << A[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    
    std::cout << "Wektor x:" << std::endl;
    std::cout << "[";
    for (size_t i = 0; i < x.size(); ++i) {
        std::cout << x[i];
        if (i < x.size() - 1) std::cout << ", ";
    }
    std::cout << "]" << std::endl;
    std::cout << std::endl;
    
    // Mnożenie gęste (dla porównania)
    std::cout << "--- Mnożenie gęste A * x ---" << std::endl;
    Vector Ax = matrixVectorMult(A, x);
    std::cout << "Wynik (A * x):" << std::endl;
    std::cout << "[";
    for (size_t i = 0; i < Ax.size(); ++i) {
        std::cout << Ax[i];
        if (i < Ax.size() - 1) std::cout << ", ";
    }
    std::cout << "]" << std::endl;
    std::cout << std::endl;
    
    // Buduj H-macierz
    int maxRank = 4;
    double epsilon = 1e-6;
    
    std::cout << "--- Kompresja do H-macierzy ---" << std::endl;
    std::cout << "Parametry: maxRank=" << maxRank << ", epsilon=" << epsilon << std::endl;
    auto H = buildHMatrix(A, maxRank, epsilon);
    std::cout << "H-macierz zbudowana." << std::endl;
    std::cout << std::endl;
    
    // Wizualizacja struktury
    std::cout << "Struktura H-macierzy (0=skompresowane, 1=podzielone):" << std::endl;
    Matrix vis = createVisualization(H);
    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < 8; ++j) {
            if (vis[i][j] == 0.0) std::cout << "█ ";
            else if (vis[i][j] == 0.5) std::cout << "▓ ";
            else std::cout << "░ ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    
    // Mnożenie H-macierz * wektor
    std::cout << "--- Mnożenie H-macierzy przez wektor ---" << std::endl;
    Vector Hx = hMatrixVectorMult(H, x);
    std::cout << "Wynik (H * x):" << std::endl;
    std::cout << "[";
    for (size_t i = 0; i < Hx.size(); ++i) {
        std::cout << Hx[i];
        if (i < Hx.size() - 1) std::cout << ", ";
    }
    std::cout << "]" << std::endl;
    std::cout << std::endl;
    
    // Oblicz błąd
    std::cout << "--- Porównanie wyników ---" << std::endl;
    std::cout << "Wynik gęsty:     [";
    for (size_t i = 0; i < Ax.size(); ++i) {
        std::cout << std::fixed << std::setprecision(2) << Ax[i];
        if (i < Ax.size() - 1) std::cout << ", ";
    }
    std::cout << "]" << std::endl;
    
    std::cout << "Wynik H-macierz: [";
    for (size_t i = 0; i < Hx.size(); ++i) {
        std::cout << std::fixed << std::setprecision(2) << Hx[i];
        if (i < Hx.size() - 1) std::cout << ", ";
    }
    std::cout << "]" << std::endl;
    
    double error = std::sqrt(vectorError(Ax, Hx));
    std::cout << std::endl;
    std::cout << "Błąd ||Ax - Hx||_2 = " << std::scientific << error << std::endl;
    
    return 0;
}
