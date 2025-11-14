#pragma once
#include <vector>
#include <memory>

using Matrix = std::vector<std::vector<double>>;

struct IMnozenie {
    virtual Matrix multiply(const Matrix& A, const Matrix& B) = 0;
    virtual ~IMnozenie() = default;
};
