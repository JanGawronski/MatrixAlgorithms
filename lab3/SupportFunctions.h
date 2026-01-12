#pragma once

#include <iostream>
#include <cstdint>
#include <vector>
#include <tuple>
#include <string>

using Matrix = std::vector<std::vector<double>>;
using Vector = std::vector<double>;

Matrix createRandomMatrix(int m, int n);
Matrix createRandomMatrix(int n);
void printSmall(const Matrix& M);

Matrix zeroMatrix(int rows, int cols);
Matrix zeroMatrix(int n);
Matrix identityMatrix(int n);
Matrix subMatrix(const Matrix &A, int row, int col, int rows, int cols);
Matrix operator+(const Matrix &A, const Matrix &B);
Matrix operator-(const Matrix &A, const Matrix &B);
Matrix operator*(const Matrix &A, const Matrix &B);
Matrix combine(const Matrix &A11, const Matrix &A12,
               const Matrix &A21, const Matrix &A22);
Matrix negate(const Matrix &A);
std::pair<bool,double> compareMatrices(const Matrix& X, const Matrix& Y, double tol);
Matrix pad(const Matrix& A, int rows, int cols);
Matrix trim(const Matrix& A, int rows, int cols);
int rows(const Matrix& M);
int cols(const Matrix& M);

std::tuple<Matrix, Vector, Matrix> svd_decomposition(const Matrix &A, int k, double epsilon = 1e-10);

double vec_dot(const Vector &a, const Vector &b);
Vector mat_vec_mul(const Matrix &M, const Vector &v);
double vec_norm(const Vector &v);
Matrix transpose_mul(const Matrix &A);
std::pair<Vector, double> power_iteration(const Matrix &A, int num_simulations = 100);
bool allclose_zero(const Matrix &M, double atol);
std::tuple<Matrix, Vector, Matrix> svd_decomposition(const Matrix &A, int k, double epsilon);


std::tuple<Matrix, Matrix, Matrix> loadImageRGB(const std::string& filename);
bool saveImageRGB(const std::string& filename, const Matrix& R, const Matrix& G, const Matrix& B);
bool saveImageGray(const std::string& filename, const Matrix& Gray);