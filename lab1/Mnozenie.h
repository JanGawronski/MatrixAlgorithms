#pragma once
#include <vector>
#include <memory>

using Matrix = std::vector<std::vector<double>>;

struct IMnozenie {
    virtual Matrix multiply(const Matrix& A, const Matrix& B) = 0;
    virtual ~IMnozenie() = default;
};

// factory tworząca implementację Binet'a (zdefiniowana w Binet.cpp)
std::unique_ptr<IMnozenie> createBinet();

// pomocnicza funkcja do generowania losowych macierzy (implementacja w Binet.cpp)
Matrix createRandomMatrix(int n);