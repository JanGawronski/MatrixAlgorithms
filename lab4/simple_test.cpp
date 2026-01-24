#include <iostream>
#include <iomanip>
#include <vector>

using Matrix = std::vector<std::vector<double>>;
using Vector = std::vector<double>;

// Proste mnożenie macierz-wektor
Vector multiply(const Matrix& A, const Vector& x) {
    int n = A.size();
    Vector result(n, 0.0);
    
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            result[i] += A[i][j] * x[j];
        }
    }
    
    return result;
}

int main() {
    std::cout << "=== Mnozenie macierzy przez wektor [1,2,3,4,5,6,7,8] ===" << std::endl;
    std::cout << std::endl;
    
    // Przykładowa macierz 8x8 (struktura z tablicy)
    Matrix A = {
        {1, 1, 0, 0, 0, 0, 0, 0},
        {1, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0}
    };
    
    // Wektor [1,2,3,4,5,6,7,8]
    Vector x = {1, 2, 3, 4, 5, 6, 7, 8};
    
    std::cout << "Macierz A:" << std::endl;
    for (int i = 0; i < 8; ++i) {
        std::cout << "  [";
        for (int j = 0; j < 8; ++j) {
            std::cout << std::setw(2) << A[i][j];
            if (j < 7) std::cout << " ";
        }
        std::cout << "]" << std::endl;
    }
    std::cout << std::endl;
    
    std::cout << "Wektor x = [1, 2, 3, 4, 5, 6, 7, 8]" << std::endl;
    std::cout << std::endl;
    
    // Mnożenie
    Vector result = multiply(A, x);
    
    std::cout << "Wynik A * x:" << std::endl;
    std::cout << "[";
    for (size_t i = 0; i < result.size(); ++i) {
        std::cout << result[i];
        if (i < result.size() - 1) std::cout << ", ";
    }
    std::cout << "]" << std::endl;
    std::cout << std::endl;
    
    std::cout << "Rozpisane:" << std::endl;
    for (size_t i = 0; i < result.size(); ++i) {
        std::cout << "  y[" << i << "] = " << result[i] << std::endl;
    }
    
    return 0;
}
