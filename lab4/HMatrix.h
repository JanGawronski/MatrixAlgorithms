#pragma once

#include <vector>
#include <memory>
#include <tuple>

using Matrix = std::vector<std::vector<double>>;
using Vector = std::vector<double>;

// H-Matrix Node structure
struct HNode {
    // For internal nodes: 4 children (top-left, top-right, bottom-left, bottom-right)
    std::vector<std::shared_ptr<HNode>> sons;
    
    // For leaf nodes: low-rank approximation M ≈ U * V^T
    int rank;
    Matrix U;  // size: rows × rank
    Matrix V;  // size: rank × cols
    
    // Matrix dimensions
    int rows;
    int cols;
    
    HNode(int r = 0, int c = 0) : rank(0), rows(r), cols(c) {}
    
    bool isLeaf() const { return sons.empty(); }
};

// Generate 3D grid topology matrix
Matrix generate3DGridMatrix(int k);

// Build H-Matrix from dense matrix using recursive compression
std::shared_ptr<HNode> buildHMatrix(const Matrix& A, int maxRank, double epsilon);

// H-Matrix operations
Vector hMatrixVectorMult(const std::shared_ptr<HNode>& H, const Vector& x);
std::shared_ptr<HNode> hMatrixAdd(const std::shared_ptr<HNode>& A, const std::shared_ptr<HNode>& B, int maxRank, double epsilon);
std::shared_ptr<HNode> hMatrixMult(const std::shared_ptr<HNode>& A, const std::shared_ptr<HNode>& B, int maxRank, double epsilon);

// Utility functions
Matrix hMatrixToDense(const std::shared_ptr<HNode>& H);
void drawHMatrix(const std::shared_ptr<HNode>& H, Matrix& visualization, int rowOffset, int colOffset);
Matrix createVisualization(const std::shared_ptr<HNode>& H);

// Helper functions from lab3
Matrix createRandomMatrix(int m, int n);
Matrix createRandomMatrix(int n);
Matrix identityMatrix(int n);
Matrix subMatrix(const Matrix& A, int startRow, int startCol, int numRows, int numCols);
Matrix combineBlocks(const Matrix& A11, const Matrix& A12, const Matrix& A21, const Matrix& A22);
Matrix matrixMultiply(const Matrix& A, const Matrix& B);
Vector matrixVectorMult(const Matrix& A, const Vector& x);
Matrix zeroMatrix(int rows, int cols);
Vector zeroVector(int size);
int rows(const Matrix& M);
int cols(const Matrix& M);
double frobeniusNorm(const Matrix& A);
double vectorError(const Vector& v1, const Vector& v2);

// Matrix operators (from lab3)
Matrix operator+(const Matrix& A, const Matrix& B);
Matrix operator-(const Matrix& A, const Matrix& B);
Matrix operator*(const Matrix& A, const Matrix& B);
Matrix negate(const Matrix& A);
std::pair<bool, double> compareMatrices(const Matrix& X, const Matrix& Y, double tol = 1e-9);
Matrix pad(const Matrix& A, int newRows, int newCols);
Matrix trim(const Matrix& A, int targetRows, int targetCols);
void printSmall(const Matrix& M);

// SVD decomposition (reuse from lab3)
std::tuple<Matrix, Vector, Matrix> svd_decomposition(const Matrix& A, int rank, double epsilon = 1e-10);

// SVD helper functions (from lab3)
double vec_dot(const Vector& a, const Vector& b);
Vector mat_vec_mul(const Matrix& M, const Vector& v);
double vec_norm(const Vector& v);
Matrix transpose_mul(const Matrix& A);
std::pair<Vector, double> power_iteration(const Matrix& A, int num_simulations = 100);
bool allclose_zero(const Matrix& M, double atol = 1e-10);

// Image I/O functions (from lab3) - for visualization
std::tuple<Matrix, Matrix, Matrix> loadImageRGB(const std::string& filename);
bool saveImageRGB(const std::string& filename, const Matrix& R, const Matrix& G, const Matrix& B);
bool saveImageGray(const std::string& filename, const Matrix& Gray);
