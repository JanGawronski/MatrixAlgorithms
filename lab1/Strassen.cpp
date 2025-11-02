#include "Strassen.h"
#include "SupportFunctions.h"

#include <algorithm>
#include <stdexcept>
#include <memory>
#include <vector>

namespace {

using Matrix = ::Matrix;

class VirtualMatrix {
    std::vector<std::vector<double>> matrix;
    int numRows;
    int numCols;
    int realRows;
    int realCols;
public:
    VirtualMatrix(const Matrix &matrix, int numRows, int numCols, int realRows, int realCols)
        :   numRows(numRows), numCols(numCols),
            realRows(realRows), realCols(realCols) {
                this->matrix = zeroMatrix(realRows, realCols);
                for (int i = 0; i < realRows; i++)
                    for (int j = 0; j < realCols; j++)
                        this->matrix[i][j] = matrix[i][j];
            }
    
        double get(int i, int j) const {
            if (i >= 0 && i < numRows && j >= 0 && j < numCols) {
                if (i < realRows && j < realCols) {
                    return matrix[i][j];
                } else {
                    return 0.0;
                }
            }
            return 0.0;
            throw std::out_of_range("MatrixView index out of range");
        }
    
        int getRows() const {
            return numRows;
        }
    
        int getCols() const {
            return numCols;
        }

        int getRealRows() const {
            return realRows;
        }

        int getRealCols() const {
            return realCols;
        }

        VirtualMatrix add(const VirtualMatrix &other, int realRows, int realCols) const {
            int Arows = std::min(this->getRealRows(), realRows);
            int Acols = std::min(this->getRealCols(), realCols);
            int Brows = std::min(other.getRealRows(), realRows);
            int Bcols = std::min(other.getRealCols(), realCols);
            Matrix result = zeroMatrix(std::max(Arows, Brows), std::max(Acols, Bcols));
            for (int i = 0; i < Arows; i++)
                for (int j = 0; j < Acols; j++)
                    result[i][j] = this->get(i, j);

            for (int i = 0; i < Brows; i++)
                for (int j = 0; j < Bcols; j++)
                    result[i][j] += other.get(i, j);

            return VirtualMatrix(result, std::max(this->getRows(), other.getRows()), std::max(this->getCols(), other.getCols()), std::max(Arows, Brows), std::max(Acols, Bcols));
        }

        VirtualMatrix substract(const VirtualMatrix &other, int realRows, int realCols) const {
            int Arows = std::min(this->getRealRows(), realRows);
            int Acols = std::min(this->getRealCols(), realCols);
            int Brows = std::min(other.getRealRows(), realRows);
            int Bcols = std::min(other.getRealCols(), realCols);
            Matrix result = zeroMatrix(std::max(Arows, Brows), std::max(Acols, Bcols));
            for (int i = 0; i < Arows; i++)
                for (int j = 0; j < Acols; j++)
                    result[i][j] = this->get(i, j);

            for (int i = 0; i < Brows; i++)
                for (int j = 0; j < Bcols; j++)
                    result[i][j] -= other.get(i, j);

            return VirtualMatrix(result, std::max(this->getRows(), other.getRows()), std::max(this->getCols(), other.getCols()), std::max(Arows, Brows), std::max(Acols, Bcols));
        }
    
        Matrix toMatrix() const {
            Matrix result = zeroMatrix(this->getRealRows(), this->getRealCols());
            for (int i = 0; i < this->getRealRows(); i++)
                for (int j = 0; j < this->getRealCols(); j++)
                    result[i][j] = this->get(i, j);
            return result;
        }

        VirtualMatrix subMatrix(int row, int col, int rows, int cols) const {
            int realSubRows = std::min(rows, this->getRealRows() - row);
            int realSubCols = std::min(cols, this->getRealCols() - col);
            Matrix sub = zeroMatrix(realSubRows, realSubCols);
            for (int i = 0; i < realSubRows; i++)
                for (int j = 0; j < realSubCols; j++)
                    sub[i][j] = this->get(row + i, col + j);
            return VirtualMatrix(sub, rows, cols, realSubRows, realSubCols);
        }

        VirtualMatrix operator*(const VirtualMatrix &other) const {
            Matrix result = zeroMatrix(this->getRows(), other.getCols());
            for (int i = 0; i < this->getCols(); i++) {
                for (int j = 0; j < other.getRows(); j++) {
                    for (int k = 0; k < this->getRows(); k++) {
                        result[k][j] += this->get(k, i) * other.get(i, j);
                    }
                }
            }
            return VirtualMatrix(result, this->getRows(), other.getCols(), this->getRealRows(), other.getRealCols());
        }
};

Matrix trim(const Matrix &A, int newRows, int newCols) {
    int Arows = static_cast<int>(A.size());
    int Acols = Arows ? static_cast<int>(A[0].size()) : 0;
    if (newRows > Arows || newCols > Acols) {
        throw std::runtime_error("New size must be smaller than original");
    }
    Matrix B = zeroMatrix(newRows, newCols);
    for (int i = 0; i < newRows; i++)
        for (int j = 0; j < newCols; j++)
            B[i][j] = A[i][j];
    return B;
}

Matrix subMatrixPadded(const Matrix &A, int row, int col, int rows, int cols, int paddedRows, int paddedCols) {
    Matrix R = zeroMatrix(paddedRows, paddedCols);
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            R[i][j] = A[row + i][col + j];
    return R;
}

VirtualMatrix multiplyRec(const VirtualMatrix &A, const VirtualMatrix &B, int realRows, int realCols) {
    int neededRows = std::min(std::min(A.getRealRows(), B.getRealCols()), realRows);
    int neededCols = std::min(std::min(A.getRealCols(), B.getRealRows()), realCols);

    if (neededRows == 1 || neededCols == 1) {
        Matrix result = zeroMatrix(neededRows, neededCols);
        for (int i = 0; i < neededRows; i++)
            for (int j = 0; j < neededCols; j++)
                for (int k = 0; k < A.getRealCols(); k++)
                    result[i][j] += A.get(i, k) * B.get(k, j);
        VirtualMatrix test = A * B;
        return VirtualMatrix(result, A.getRows(), B.getCols(), neededRows, neededCols);
    }

    if (neededRows == 0 || neededCols == 0) {
        return VirtualMatrix(zeroMatrix(0, 0), A.getRows(), B.getCols(), neededRows, neededCols);
    }

    int C11rows = (A.getRows() + 1) / 2;
    int C11cols = (A.getCols() + 1) / 2;
    int C12rows = C11rows;
    int C12cols = B.getCols() - C11cols;
    int C21rows = A.getRows() - C11rows;
    int C21cols = C11cols;
    int C22rows = C21rows;
    int C22cols = C12cols;

    int C11realrows = std::min(C11rows, neededRows);
    int C11realcols = std::min(C11cols, neededCols);

    int C12realrows = C11realrows;
    int C12realcols = std::max(neededCols - C11cols, 0);

    int C21realrows = std::max(neededRows - C11rows, 0);
    int C21realcols = C11realcols;

    int C22realrows = C21realrows;
    int C22realcols = C12realcols;

    int P1rows = C11realrows;
    int P1cols = C11realcols;

    int P2rows = C21realrows;
    int P2cols = C21realcols;

    int P3rows = C12realrows;
    int P3cols = C12realcols;

    int P4rows = C11realrows;
    int P4cols = C11realcols;

    int P5rows = C11realrows;
    int P5cols = C11realcols;

    int P6rows = C22realrows;
    int P6cols = C22realcols;

    int P7rows = C11realrows;
    int P7cols = C11realcols;

    VirtualMatrix C11 = VirtualMatrix(zeroMatrix(P1rows, P1cols), C11rows, C11cols, P1rows, P1cols);
    VirtualMatrix C12 = VirtualMatrix(zeroMatrix(P3rows, P3cols), C12rows, C12cols, P3rows, P3cols);
    VirtualMatrix C21 = VirtualMatrix(zeroMatrix(P2rows, P2cols), C21rows, C21cols, P2rows, P2cols);
    VirtualMatrix C22 = VirtualMatrix(zeroMatrix(P6rows, P6cols), C22rows, C22cols, P6rows, P6cols);

    VirtualMatrix A11 = A.subMatrix(0, 0, C11rows, C11cols);
    VirtualMatrix A12 = A.subMatrix(0, C11cols, C12rows, C12cols);
    VirtualMatrix A21 = A.subMatrix(C11rows, 0, C21rows, C21cols);
    VirtualMatrix A22 = A.subMatrix(C11rows, C11cols, C22rows, C22cols);

    VirtualMatrix B11 = B.subMatrix(0, 0, C11cols, C11rows);
    VirtualMatrix B12 = B.subMatrix(0, C11rows, C12cols, C11cols);
    VirtualMatrix B21 = B.subMatrix(C11cols, 0, C21cols, C11rows);
    VirtualMatrix B22 = B.subMatrix(C11cols, C11rows, C22cols, C22rows);

    if (P1rows > 0 && P1cols > 0) {
        VirtualMatrix P1 = multiplyRec(A11.add(A22, P1rows, P1cols), B11.add(B22, P1rows, P1cols), P1rows, P1cols);
        C11 = C11.add(P1, C11rows, C11cols);
        C22 = C22.add(P1, C22rows, C22cols);
    }

    if (P2rows > 0 && P2cols > 0) {
        VirtualMatrix P2 = multiplyRec(A21.add(A22, P2rows, P2cols), B11, P2rows, P2cols);
        C21 = C21.add(P2, C21rows, C21cols);
        C22 = C22.substract(P2, C22rows, C22cols);
    }

    if (P3rows > 0 && P3cols > 0) {
        VirtualMatrix P3 = multiplyRec(A11, B12.substract(B22, P3rows, P3cols), P3rows, P3cols);
        C12 = C12.add(P3, C12rows, C12cols);
        C22 = C22.add(P3, C22rows, C22cols);
    }

    if (P4rows > 0 && P4cols > 0) {
        VirtualMatrix P4 = multiplyRec(A22, B21.substract(B11, P4rows, P4cols), P4rows, P4cols);
        C11 = C11.add(P4, C21rows, C21cols);
        C21 = C21.add(P4, C22rows, C22cols);
    }

    if (P5rows > 0 && P5cols > 0) {
        VirtualMatrix P5 = multiplyRec(A11.add(A12, P5rows, P5cols), B22, P5rows, P5cols);
        C11 = C11.substract(P5, C11rows, C11cols);
        C12 = C12.add(P5, C12rows, C12cols);
    }

    if (P6rows > 0 && P6cols > 0) {
        VirtualMatrix P6 = multiplyRec(A21.substract(A11, P6rows, P6cols), B11.add(B12, P6rows, P6cols), P6rows, P6cols);
        C22 = C22.add(P6, C22rows, C22cols);
    }

    if (P7rows > 0 && P7cols > 0) {
        VirtualMatrix P7 = multiplyRec(A12.substract(A22, P7rows, P7cols), B21.add(B22, P7rows, P7cols), P7rows, P7cols);
        C11 = C11.add(P7, C11rows, C11cols);
    }

    Matrix result = zeroMatrix(neededRows, neededCols);
    for (int i = 0; i < C11realrows; i++)
        for (int j = 0; j < C11realcols; j++)
            result[i][j] = C11.get(i, j);

    for (int i = 0; i < C12realrows; i++)
        for (int j = 0; j < C12realcols; j++)
            result[i][j + C11cols] = C12.get(i, j);

    for (int i = 0; i < C21realrows; i++)
        for (int j = 0; j < C21realcols; j++)
            result[i + C11rows][j] = C21.get(i, j);

    for (int i = 0; i < C22realrows; i++)
        for (int j = 0; j < C22realcols; j++)
            result[i + C11rows][j + C11cols] = C22.get(i, j);

    VirtualMatrix test = A * B;

    return VirtualMatrix(result, A.getRows(), B.getCols(), neededRows, neededCols);

}

}

class StrassenImpl : public IMnozenie {
public:
    Matrix multiply(const Matrix &A, const Matrix &B) override {
        VirtualMatrix vA(A, static_cast<int>(A.size()), static_cast<int>(A[0].size()), static_cast<int>(A.size()), static_cast<int>(A[0].size()));
        VirtualMatrix vB(B, static_cast<int>(B.size()), static_cast<int>(B[0].size()), static_cast<int>(B.size()), static_cast<int>(B[0].size()));
        return multiplyRec(vA, vB, static_cast<int>(A.size()), static_cast<int>(B[0].size())).toMatrix();
    }
};

std::unique_ptr<IMnozenie> createStrassen() {
    return std::make_unique<StrassenImpl>();
}

