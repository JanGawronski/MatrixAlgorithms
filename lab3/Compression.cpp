#include "Compression.h"
#include <algorithm>

TreeNode* createTree(const Matrix &A, int rank, double epsilon) {
    auto [U, D, V] = svd_decomposition(A, rank);
   
    if (static_cast<int>(D.size()) < rank || D[rank - 1] < epsilon) {
        if (allclose_zero(A, 1e-10)) {
            TreeNode* node = new TreeNode();
            node->singularValues = zeroMatrix(1, rank)[0];
            node->U = zeroMatrix(rows(A), rank);
            node->V = zeroMatrix(rank, cols(A));
            return node;
        }
        TreeNode* node = new TreeNode();
        node->singularValues = D;
        node->U = subMatrix(U, 0, 0, rows(U), static_cast<int>(D.size()));
        node->V = subMatrix(V, 0, 0, static_cast<int>(D.size()), cols(V));
        return node;
    }
    TreeNode* node = new TreeNode();
    int rmid = rows(A) / 2;
    int cmid = cols(A) / 2;
    node->topLeft = createTree(subMatrix(A, 0, 0, rmid, cmid), rank, epsilon);
    node->topRight = createTree(subMatrix(A, 0, cmid, rmid, cols(A) - cmid), rank, epsilon);
    node->bottomLeft = createTree(subMatrix(A, rmid, 0, rows(A) - rmid, cmid), rank, epsilon);
    node->bottomRight = createTree(subMatrix(A, rmid, cmid, rows(A) - rmid, cols(A) - cmid), rank, epsilon);    
    return node;
}

Matrix reconstructFromTree(TreeNode* node) {
    if (node->topLeft == nullptr && node->topRight == nullptr &&
        node->bottomLeft == nullptr && node->bottomRight == nullptr) {
        int r = static_cast<int>(node->singularValues.size());
        Matrix S = zeroMatrix(r, r);
        for (int i = 0; i < r; ++i) S[i][i] = node->singularValues[i];
        return node->U * S * node->V;
    }
    Matrix A11 = reconstructFromTree(node->topLeft);
    Matrix A12 = reconstructFromTree(node->topRight);
    Matrix A21 = reconstructFromTree(node->bottomLeft);
    Matrix A22 = reconstructFromTree(node->bottomRight);
    return combine(A11, A12, A21, A22);
}

void drawCompression(TreeNode* node, Matrix &matrix, int x0, int x1, int y0, int y1) {
    if (node == nullptr) return;
     if (node->topLeft == nullptr && node->topRight == nullptr &&
        node->bottomLeft == nullptr && node->bottomRight == nullptr) {
        int rank = static_cast<int>(node->singularValues.size());
        if (rank <= 0) return;
        int maxRows = rows(matrix);
        int maxCols = cols(matrix);
        int x_end = std::clamp(x1, 0, maxRows);
        int y_end = std::clamp(y1, 0, maxCols);
        int y_rank_end = std::min(y_end, y0 + rank);
        int x_rank_end = std::min(x_end, x0 + rank);
        for (int i = std::clamp(x0, 0, maxRows); i < x_end; ++i) {
            for (int j = std::clamp(y0, 0, maxCols); j < y_rank_end; ++j)
                matrix[i][j] = 0.0;
        }
        for (int i = std::clamp(x0, 0, maxRows); i < x_rank_end; ++i) {
            for (int j = std::clamp(y0, 0, maxCols); j < y_end; ++j)
                matrix[i][j] = 0.0;
        }
        return;
    }

    int x_mid = (x0 + x1) / 2;
    int y_mid = (y0 + y1) / 2;
    drawCompression(node->topLeft, matrix, x0, x_mid, y0, y_mid);
    drawCompression(node->topRight, matrix, x0, x_mid, y_mid, y1);
    drawCompression(node->bottomLeft, matrix, x_mid, x1, y0, y_mid);
    drawCompression(node->bottomRight, matrix, x_mid, x1, y_mid, y1);
}

Matrix drawCompression(TreeNode* node, int width, int height) {
    Matrix matrix = zeroMatrix(width, height);
    for (int i = 0; i < width; ++i)
        for (int j = 0; j < height; ++j)
            matrix[i][j] = 1.0;
    drawCompression(node, matrix, 0, rows(matrix), 0, cols(matrix));
    return matrix;
}