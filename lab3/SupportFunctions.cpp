#include "SupportFunctions.h"
#include <random>
#include <iomanip>
#include <vector>
#include "stb_image.h"
#include "stb_image_write.h"
#include <cstring>

Matrix createRandomMatrix(int m, int n) {
    Matrix M(m, std::vector<double>(n));
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

void printSmall(const Matrix& M) {
    for (int i = 0; i < rows(M); ++i) {
        for (int j = 0; j < cols(M); ++j)
            std::cout << std::setprecision(6) << std::setw(12) << M[i][j];
        std::cout << '\n';
    }
}

Matrix zeroMatrix(int rows, int cols) {
    return Matrix(rows, std::vector<double>(cols, 0.0));
}

Matrix zeroMatrix(int n) {
    return zeroMatrix(n, n);
}

// A[row:row+rows-1][col:col+cols-1]
Matrix subMatrix(const Matrix &A, int row, int col, int rows, int cols) {
    Matrix R = zeroMatrix(rows, cols);
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            R[i][j] = A[row + i][col + j];
    return R;
}

Matrix combine(const Matrix &A11, const Matrix &A12,
               const Matrix &A21, const Matrix &A22) {
    Matrix A = zeroMatrix(rows(A11) + rows(A21), cols(A11) + cols(A12));
    for (int i = 0; i < rows(A11); i++)
        for (int j = 0; j < cols(A11); j++)
            A[i][j] = A11[i][j];

    for (int i = 0; i < rows(A12); i++)
        for (int j = 0; j < cols(A12); j++)
            A[i][j + cols(A11)] = A12[i][j];

    for (int i = 0; i < rows(A21); i++)
        for (int j = 0; j < cols(A11); j++)
            A[i + rows(A11)][j] = A21[i][j];

    for (int i = 0; i < rows(A22); i++)
        for (int j = 0; j < cols(A12); j++)
            A[i + rows(A11)][j + cols(A11)] = A22[i][j];

    return A;
}

Matrix pad(const Matrix& A, int newRows, int newCols) {
    Matrix P = zeroMatrix(newRows, newCols);
    for (int i = 0; i < rows(A); ++i)
        for (int j = 0; j < cols(A); ++j)
            P[i][j] = A[i][j];
    return P;
}

Matrix trim(const Matrix& A, int rows, int cols) {
    Matrix T = zeroMatrix(rows, cols);
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            T[i][j] = A[i][j];
    return T;
}

int rows(const Matrix& M) {
    return static_cast<int>(M.size());
}

int cols(const Matrix& M) {
    return M.size() ? static_cast<int>(M[0].size()) : 0;
}

Matrix operator+(const Matrix &A, const Matrix &B) {
    Matrix R = zeroMatrix(std::max(rows(A), rows(B)), std::max(cols(A), cols(B)));
    for (int i = 0; i < rows(A); ++i)
        for (int j = 0; j < cols(A); ++j)
            R[i][j] = A[i][j];

    for (int i = 0; i < rows(B); ++i) {
        for (int j = 0; j < cols(B); ++j)
            R[i][j] += B[i][j];

    }

    return R;
}

Matrix operator-(const Matrix &A, const Matrix &B) {
    Matrix R = zeroMatrix(std::max(rows(A), rows(B)), std::max(cols(A), cols(B)));
    for (int i = 0; i < rows(A); ++i)
        for (int j = 0; j < cols(A); ++j)
            R[i][j] = A[i][j];

    for (int i = 0; i < rows(B); ++i) {
        for (int j = 0; j < cols(B); ++j)
            R[i][j] -= B[i][j];
    }

    return R;
}

Matrix operator*(const Matrix &A, const Matrix &B) {

    Matrix C = zeroMatrix(rows(A), cols(B));
    for (int i = 0; i < rows(A); ++i) {
        for (int k = 0; k < cols(A); ++k) {
            double aik = A[i][k];
            for (int j = 0; j < cols(B); ++j) {
                double prod = aik * B[k][j];
                C[i][j] += prod;
            }
        }
    }

    return C;
}

Matrix negate(const Matrix &A) {
    Matrix R = zeroMatrix(rows(A), cols(A));
    for (int i = 0; i < rows(A); ++i)
        for (int j = 0; j < cols(A); ++j)
            R[i][j] = -A[i][j];
    return R;
}

std::pair<bool,double> compareMatrices(const Matrix& X, const Matrix& Y, double tol) {
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


double vec_dot(const Vector &a, const Vector &b) {
    return std::inner_product(a.begin(), a.end(), b.begin(), 0.0);
}

Vector mat_vec_mul(const Matrix &M, const Vector &v) {
    size_t cols = v.size();
    Vector r(rows(M), 0.0);
    for (int i = 0; i < rows(M); ++i) {
        double s = 0.0;
        for (size_t j = 0; j < cols; ++j) s += M[i][j] * v[j];
        r[i] = s;
    }
    return r;
}

double vec_norm(const Vector &v) {
    return std::sqrt(vec_dot(v, v));
}

Matrix transpose_mul(const Matrix &A) {
    size_t m = A.size(), n = A.empty() ? 0 : A[0].size();
    Matrix B(n, Vector(n, 0.0));
    for (size_t i = 0; i < m; ++i)
        for (size_t j = 0; j < n; ++j)
            for (size_t k = 0; k < n; ++k)
                B[j][k] += A[i][j] * A[i][k];
    return B;
}

std::pair<Vector, double> power_iteration(const Matrix &A, int num_simulations) {
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

bool allclose_zero(const Matrix &M, double atol) {
    for (const auto &row : M)
        for (double v : row)
            if (std::abs(v) > atol) return false;
    return true;
}

std::tuple<Matrix, Vector, Matrix> svd_decomposition(const Matrix &A, int r, double epsilon) {
    size_t m = A.size(), n = A.empty() ? 0 : A[0].size();
    Matrix U(m, Vector(r, 0.0));
    Matrix V(r, Vector(n, 0.0));
    Vector S;
    Matrix B = transpose_mul(A);
    int num_valid = 0;

    for (int i = 0; i < r; ++i) {
        auto [v, sigma_squared] = power_iteration(B);
        if (sigma_squared < epsilon * epsilon) break;
        double sigma = std::sqrt(sigma_squared);
        S.push_back(sigma);
        if ((size_t)i < V.size()) V[i] = v;
        Vector u = mat_vec_mul(A, v);
        if (sigma != 0.0) for (double &x : u) x /= sigma;
        for (size_t r = 0; r < m; ++r) U[r][num_valid] = u[r];
        for (size_t p = 0; p < n; ++p)
            for (size_t q = 0; q < n; ++q)
                B[p][q] -= sigma_squared * v[p] * v[q];
        ++num_valid;
        if (allclose_zero(B, epsilon)) break;
    }

    Matrix U_trim(m, Vector(num_valid));
    Matrix V_trim(num_valid, Vector(n));
    for (size_t r = 0; r < m; ++r)
        for (int c = 0; c < num_valid; ++c)
            U_trim[r][c] = U[r][c];
    for (int r = 0; r < num_valid; ++r)
        for (size_t c = 0; c < n; ++c)
            V_trim[r][c] = V[r][c];

    return {U_trim, S, V_trim};
}

std::tuple<Matrix, Matrix, Matrix> loadImageRGB(const std::string& filename) {
    int w=0, h=0, c=0;
    unsigned char *data = stbi_load(filename.c_str(), &w, &h, &c, 3); // force RGB
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
    unsigned char *buf = new unsigned char[bufSize];
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

bool save
