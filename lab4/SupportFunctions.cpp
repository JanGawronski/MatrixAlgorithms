#include "HMatrix.h"
#include <cmath>
#include <random>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <limits>

// STB Image libraries for visualization
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

// ============================================================================
// Matrix Creation Functions (from lab3)
// ============================================================================

Matrix createRandomMatrix(int m, int n) {
    Matrix M(m, Vector(n));
    std::mt19937_64 rng(std::random_device{}());
    std::uniform_real_distribution<double> dist(1e-8 + 1e-16, 1.0 - 1e-16);

    for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j)
            M[i][j] = dist(rng);
    return M;
}

Matrix createRandomMatrix(int n) {
    return createRandomMatrix(n, n);
}

Matrix identityMatrix(int n) {
    Matrix I = zeroMatrix(n, n);
    for (int i = 0; i < n; ++i)
        I[i][i] = 1.0;
    return I;
}

// ============================================================================
// Matrix Helper Functions
// ============================================================================

Matrix subMatrix(const Matrix& A, int startRow, int startCol, int numRows, int numCols) {
    Matrix result(numRows, Vector(numCols));
    for (int i = 0; i < numRows; ++i) {
        for (int j = 0; j < numCols; ++j) {
            result[i][j] = A[startRow + i][startCol + j];
        }
    }
    return result;
}

Matrix combineBlocks(const Matrix& A11, const Matrix& A12, 
                     const Matrix& A21, const Matrix& A22) {
    int rows1 = rows(A11);
    int rows2 = rows(A21);
    int cols1 = cols(A11);
    int cols2 = cols(A12);
    
    Matrix result(rows1 + rows2, Vector(cols1 + cols2));
    
    // Top-left
    for (int i = 0; i < rows1; ++i) {
        for (int j = 0; j < cols1; ++j) {
            result[i][j] = A11[i][j];
        }
    }
    
    // Top-right
    for (int i = 0; i < rows1; ++i) {
        for (int j = 0; j < cols2; ++j) {
            result[i][cols1 + j] = A12[i][j];
        }
    }
    
    // Bottom-left
    for (int i = 0; i < rows2; ++i) {
        for (int j = 0; j < cols1; ++j) {
            result[rows1 + i][j] = A21[i][j];
        }
    }
    
    // Bottom-right
    for (int i = 0; i < rows2; ++i) {
        for (int j = 0; j < cols2; ++j) {
            result[rows1 + i][cols1 + j] = A22[i][j];
        }
    }
    
    return result;
}

Matrix matrixMultiply(const Matrix& A, const Matrix& B) {
    int m = rows(A);
    int n = cols(B);
    int k = cols(A);
    
    if (k != rows(B)) {
        throw std::runtime_error("Matrix dimensions don't match for multiplication");
    }
    
    Matrix result = zeroMatrix(m, n);
    
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            double sum = 0.0;
            for (int p = 0; p < k; ++p) {
                sum += A[i][p] * B[p][j];
            }
            result[i][j] = sum;
        }
    }
    
    return result;
}

Vector matrixVectorMult(const Matrix& A, const Vector& x) {
    int m = rows(A);
    int n = cols(A);
    
    if (n != static_cast<int>(x.size())) {
        throw std::runtime_error("Matrix-vector dimensions don't match");
    }
    
    Vector result(m, 0.0);
    
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            result[i] += A[i][j] * x[j];
        }
    }
    
    return result;
}

Matrix zeroMatrix(int rows, int cols) {
    return Matrix(rows, Vector(cols, 0.0));
}

Vector zeroVector(int size) {
    return Vector(size, 0.0);
}

int rows(const Matrix& M) {
    return static_cast<int>(M.size());
}

int cols(const Matrix& M) {
    return M.empty() ? 0 : static_cast<int>(M[0].size());
}

double frobeniusNorm(const Matrix& A) {
    double sum = 0.0;
    for (const auto& row : A) {
        for (double val : row) {
            sum += val * val;
        }
    }
    return std::sqrt(sum);
}

double vectorError(const Vector& v1, const Vector& v2) {
    if (v1.size() != v2.size()) {
        throw std::runtime_error("Vector sizes don't match");
    }
    
    double sum = 0.0;
    for (size_t i = 0; i < v1.size(); ++i) {
        double diff = v1[i] - v2[i];
        sum += diff * diff;
    }
    return sum;
}

// ============================================================================
// SVD Helper Functions (from lab3)
// ============================================================================

double vec_dot(const Vector& a, const Vector& b) {
    double result = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        result += a[i] * b[i];
    }
    return result;
}

Vector mat_vec_mul(const Matrix& M, const Vector& v) {
    return matrixVectorMult(M, v);
}

double vec_norm(const Vector& v) {
    return std::sqrt(vec_dot(v, v));
}

Matrix transpose_mul(const Matrix& A) {
    int m = rows(A);
    int n = cols(A);
    Matrix result = zeroMatrix(n, n);
    
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < m; ++k) {
                result[i][j] += A[k][i] * A[k][j];
            }
        }
    }
    
    return result;
}

bool allclose_zero(const Matrix& M, double atol) {
    for (const auto& row : M) {
        for (double val : row) {
            if (std::abs(val) > atol) return false;
        }
    }
    return true;
}

std::pair<Vector, double> power_iteration(const Matrix& A, int num_simulations) {
    size_t n = A.size();
    std::mt19937_64 gen(std::random_device{}());
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    Vector b_k(n);
    for (size_t i = 0; i < n; ++i) b_k[i] = dist(gen);

    for (int it = 0; it < num_simulations; ++it) {
        Vector b_k1 = mat_vec_mul(A, b_k);
        double norm = vec_norm(b_k1);
        if (norm == 0.0) break;
        for (size_t i = 0; i < n; ++i) b_k[i] = b_k1[i] / norm;
    }
    
    Vector Ab = mat_vec_mul(A, b_k);
    double denom = vec_dot(b_k, b_k);
    double eigenvalue = denom == 0.0 ? 0.0 : vec_dot(b_k, Ab) / denom;
    return {b_k, eigenvalue};
}

std::tuple<Matrix, Vector, Matrix> svd_decomposition(const Matrix& A, int r, double epsilon) {
    size_t m = A.size();
    size_t n = A.empty() ? 0 : A[0].size();
    
    Matrix U(m, Vector(r, 0.0));
    Matrix V(r, Vector(n, 0.0));
    Vector S;
    Matrix B = transpose_mul(A);
    int num_valid = 0;

    for (int i = 0; i < r; ++i) {
        auto [v, sigma_squared] = power_iteration(B, 100);
        if (sigma_squared < epsilon * epsilon) break;
        double sigma = std::sqrt(sigma_squared);
        S.push_back(sigma);
        if ((size_t)i < V.size()) V[i] = v;
        Vector u = mat_vec_mul(A, v);
        if (sigma != 0.0) for (double& x : u) x /= sigma;
        for (size_t r_idx = 0; r_idx < m; ++r_idx) U[r_idx][num_valid] = u[r_idx];
        for (size_t p = 0; p < n; ++p)
            for (size_t q = 0; q < n; ++q)
                B[p][q] -= sigma_squared * v[p] * v[q];
        ++num_valid;
        if (allclose_zero(B, epsilon)) break;
    }

    Matrix U_trim(m, Vector(num_valid));
    Matrix V_trim(num_valid, Vector(n));
    for (size_t r_idx = 0; r_idx < m; ++r_idx)
        for (int c = 0; c < num_valid; ++c)
            U_trim[r_idx][c] = U[r_idx][c];
    for (int r_idx = 0; r_idx < num_valid; ++r_idx)
        for (size_t c = 0; c < n; ++c)
            V_trim[r_idx][c] = V[r_idx][c];

    return {U_trim, S, V_trim};
}

// ============================================================================
// Matrix Operators (from lab3)
// ============================================================================

Matrix operator+(const Matrix& A, const Matrix& B) {
    if (rows(A) != rows(B) || cols(A) != cols(B)) {
        throw std::runtime_error("Matrix dimensions don't match for addition");
    }
    Matrix R = zeroMatrix(rows(A), cols(A));
    for (int i = 0; i < rows(A); ++i)
        for (int j = 0; j < cols(A); ++j)
            R[i][j] = A[i][j] + B[i][j];
    return R;
}

Matrix operator-(const Matrix& A, const Matrix& B) {
    if (rows(A) != rows(B) || cols(A) != cols(B)) {
        throw std::runtime_error("Matrix dimensions don't match for subtraction");
    }
    Matrix R = zeroMatrix(rows(A), cols(A));
    for (int i = 0; i < rows(A); ++i)
        for (int j = 0; j < cols(A); ++j)
            R[i][j] = A[i][j] - B[i][j];
    return R;
}

Matrix operator*(const Matrix& A, const Matrix& B) {
    return matrixMultiply(A, B);
}

Matrix negate(const Matrix& A) {
    Matrix R = zeroMatrix(rows(A), cols(A));
    for (int i = 0; i < rows(A); ++i)
        for (int j = 0; j < cols(A); ++j)
            R[i][j] = -A[i][j];
    return R;
}

std::pair<bool, double> compareMatrices(const Matrix& X, const Matrix& Y, double tol) {
    if (rows(X) != rows(Y)) return {false, std::numeric_limits<double>::infinity()};
    if (rows(X) == 0) return {true, 0.0};
    if (cols(X) != cols(Y)) return {false, std::numeric_limits<double>::infinity()};
    double maxDiff = 0.0;
    for (int i = 0; i < rows(X); ++i) {
        for (int j = 0; j < cols(X); ++j) {
            double d = std::abs(X[i][j] - Y[i][j]);
            if (d > maxDiff) maxDiff = d;
            if (maxDiff > tol) return {false, maxDiff};
        }
    }
    return {true, maxDiff};
}

Matrix pad(const Matrix& A, int newRows, int newCols) {
    Matrix P = zeroMatrix(newRows, newCols);
    for (int i = 0; i < rows(A); ++i)
        for (int j = 0; j < cols(A); ++j)
            P[i][j] = A[i][j];
    return P;
}

Matrix trim(const Matrix& A, int targetRows, int targetCols) {
    Matrix T = zeroMatrix(targetRows, targetCols);
    for (int i = 0; i < targetRows; ++i)
        for (int j = 0; j < targetCols; ++j)
            T[i][j] = A[i][j];
    return T;
}

void printSmall(const Matrix& M) {
    for (int i = 0; i < rows(M); ++i) {
        for (int j = 0; j < cols(M); ++j)
            std::cout << std::setw(12) << std::setprecision(6) << M[i][j];
        std::cout << '\n';
    }
}

// ============================================================================
// Image I/O Functions (from lab3)
// ============================================================================

std::tuple<Matrix, Matrix, Matrix> loadImageRGB(const std::string& filename) {
    int w = 0, h = 0, c = 0;
    unsigned char* data = stbi_load(filename.c_str(), &w, &h, &c, 3); // force RGB
    if (!data) return {Matrix(), Matrix(), Matrix()};

    Matrix R(h, Vector(w, 0.0));
    Matrix G(h, Vector(w, 0.0));
    Matrix B(h, Vector(w, 0.0));
    size_t idx = 0;
    for (int y = 0; y < h; ++y) {
        for (int x = 0; x < w; ++x) {
            R[y][x] = static_cast<double>(data[idx++]) / 255.0;
            G[y][x] = static_cast<double>(data[idx++]) / 255.0;
            B[y][x] = static_cast<double>(data[idx++]) / 255.0;
        }
    }
    stbi_image_free(data);
    return {R, G, B};
}

bool saveImageRGB(const std::string& filename, const Matrix& R, const Matrix& G, const Matrix& B) {
    if (rows(R) == 0 || rows(R) != rows(G) || rows(R) != rows(B)) return false;
    if (cols(R) == 0 || cols(R) != cols(G) || cols(R) != cols(B)) return false;

    int h = rows(R);
    int w = cols(R);
    size_t bufSize = static_cast<size_t>(w) * static_cast<size_t>(h) * 3;
    unsigned char* buf = new unsigned char[bufSize];
    size_t idx = 0;
    for (int y = 0; y < h; ++y) {
        for (int x = 0; x < w; ++x) {
            double rv = R[y][x];
            double gv = G[y][x];
            double bv = B[y][x];
            int ri = static_cast<int>(std::round(std::clamp(rv, 0.0, 1.0) * 255.0));
            int gi = static_cast<int>(std::round(std::clamp(gv, 0.0, 1.0) * 255.0));
            int bi = static_cast<int>(std::round(std::clamp(bv, 0.0, 1.0) * 255.0));
            buf[idx++] = static_cast<unsigned char>(ri);
            buf[idx++] = static_cast<unsigned char>(gi);
            buf[idx++] = static_cast<unsigned char>(bi);
        }
    }

    int ok = stbi_write_png(filename.c_str(), w, h, 3, buf, w * 3);
    delete[] buf;
    return ok != 0;
}

bool saveImageGray(const std::string& filename, const Matrix& Gray) {
    if (rows(Gray) == 0 || cols(Gray) == 0) return false;

    int h = rows(Gray);
    int w = cols(Gray);
    size_t bufSize = static_cast<size_t>(w) * static_cast<size_t>(h);
    unsigned char* buf = new unsigned char[bufSize];
    size_t idx = 0;
    for (int y = 0; y < h; ++y) {
        for (int x = 0; x < w; ++x) {
            double gv = Gray[y][x];
            int gi = static_cast<int>(std::round(std::clamp(gv, 0.0, 1.0) * 255.0));
            buf[idx++] = static_cast<unsigned char>(gi);
        }
    }

    int ok = stbi_write_png(filename.c_str(), w, h, 1, buf, w);
    delete[] buf;
    return ok != 0;
}

