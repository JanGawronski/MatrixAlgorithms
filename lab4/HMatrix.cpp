#include "HMatrix.h"
#include <cmath>
#include <algorithm>
#include <random>
#include <iostream>

// ============================================================================
// 3D Grid Matrix Generation
// ============================================================================

Matrix generate3DGridMatrix(int k) {
    int n = 1 << (3 * k);  // 2^(3k) = (2^k)^3
    int gridSize = 1 << k; // 2^k per dimension
    
    Matrix A = zeroMatrix(n, n);
    std::mt19937 gen(42);  // Fixed seed for reproducibility
    std::uniform_real_distribution<> dis(1.0, 10.0);
    
    // Convert 1D index to 3D coordinates
    auto idx_to_3d = [gridSize](int idx) -> std::tuple<int, int, int> {
        int x = idx / (gridSize * gridSize);
        int y = (idx / gridSize) % gridSize;
        int z = idx % gridSize;
        return {x, y, z};
    };
    
    // Convert 3D coordinates to 1D index
    auto _3d_to_idx = [gridSize](int x, int y, int z) -> int {
        return x * gridSize * gridSize + y * gridSize + z;
    };
    
    // For each vertex in the grid
    for (int i = 0; i < n; ++i) {
        auto [x, y, z] = idx_to_3d(i);
        
        // Add connections to neighboring vertices (6-connectivity)
        std::vector<std::tuple<int, int, int>> neighbors = {
            {x-1, y, z}, {x+1, y, z},  // x-axis neighbors
            {x, y-1, z}, {x, y+1, z},  // y-axis neighbors
            {x, y, z-1}, {x, y, z+1}   // z-axis neighbors
        };
        
        for (auto [nx, ny, nz] : neighbors) {
            if (nx >= 0 && nx < gridSize && 
                ny >= 0 && ny < gridSize && 
                nz >= 0 && nz < gridSize) {
                int j = _3d_to_idx(nx, ny, nz);
                A[i][j] = dis(gen);  // Random weight
            }
        }
        
        // Add diagonal dominance for stability
        double rowSum = 0.0;
        for (int j = 0; j < n; ++j) {
            if (i != j) rowSum += std::abs(A[i][j]);
        }
        A[i][i] = rowSum + dis(gen);
    }
    
    return A;
}

// ============================================================================
// H-Matrix Construction
// ============================================================================

std::shared_ptr<HNode> buildHMatrix(const Matrix& A, int maxRank, double epsilon) {
    int m = rows(A);
    int n = cols(A);
    
    if (m == 0 || n == 0) {
        return std::make_shared<HNode>(m, n);
    }
    
    auto node = std::make_shared<HNode>(m, n);
    
    // Try low-rank approximation
    auto [U, S, V] = svd_decomposition(A, maxRank, epsilon);
    
    int actualRank = static_cast<int>(S.size());
    
    // Check if low-rank approximation is sufficient
    bool isLowRank = (actualRank < maxRank) || 
                     (actualRank > 0 && S[actualRank - 1] < epsilon);
    
    if (isLowRank || m <= 2 || n <= 2) {
        // Create leaf node with low-rank approximation
        node->rank = actualRank;
        node->U = U;
        node->V = V;
        return node;
    }
    
    // Subdivide into 4 blocks
    int midRow = m / 2;
    int midCol = n / 2;
    
    node->sons.resize(4);
    node->sons[0] = buildHMatrix(subMatrix(A, 0, 0, midRow, midCol), maxRank, epsilon);
    node->sons[1] = buildHMatrix(subMatrix(A, 0, midCol, midRow, n - midCol), maxRank, epsilon);
    node->sons[2] = buildHMatrix(subMatrix(A, midRow, 0, m - midRow, midCol), maxRank, epsilon);
    node->sons[3] = buildHMatrix(subMatrix(A, midRow, midCol, m - midRow, n - midCol), maxRank, epsilon);
    
    return node;
}

// ============================================================================
// H-Matrix Vector Multiplication (Slide 20)
// ============================================================================

Vector hMatrixVectorMult(const std::shared_ptr<HNode>& H, const Vector& x) {
    if (!H || H->rows == 0) {
        return zeroVector(H ? H->rows : 0);
    }
    
    // Leaf case: Y = U * (V * x)
    if (H->isLeaf()) {
        if (H->rank == 0) {
            return zeroVector(H->rows);
        }
        
        // First: temp = V * x (rank × cols) * (cols × 1) = (rank × 1)
        Vector temp = matrixVectorMult(H->V, x);
        
        // Then: Y = U * temp (rows × rank) * (rank × 1) = (rows × 1)
        return matrixVectorMult(H->U, temp);
    }
    
    // Internal node case: split vector and recurse
    int splitPoint = H->sons[0]->cols;
    
    Vector x1(x.begin(), x.begin() + splitPoint);
    Vector x2(x.begin() + splitPoint, x.end());
    
    Vector y1_1 = hMatrixVectorMult(H->sons[0], x1);
    Vector y1_2 = hMatrixVectorMult(H->sons[1], x2);
    Vector y2_1 = hMatrixVectorMult(H->sons[2], x1);
    Vector y2_2 = hMatrixVectorMult(H->sons[3], x2);
    
    // Combine: [y1_1 + y1_2; y2_1 + y2_2]
    Vector result;
    result.reserve(H->rows);
    
    for (size_t i = 0; i < y1_1.size(); ++i) {
        result.push_back(y1_1[i] + y1_2[i]);
    }
    for (size_t i = 0; i < y2_1.size(); ++i) {
        result.push_back(y2_1[i] + y2_2[i]);
    }
    
    return result;
}

// ============================================================================
// H-Matrix Addition (Slides 21-22)
// ============================================================================

std::shared_ptr<HNode> hMatrixAdd(const std::shared_ptr<HNode>& A, 
                                   const std::shared_ptr<HNode>& B, 
                                   int maxRank, double epsilon) {
    if (!A || !B || A->rows != B->rows || A->cols != B->cols) {
        throw std::runtime_error("Invalid matrix dimensions for addition");
    }
    
    auto result = std::make_shared<HNode>(A->rows, A->cols);
    
    // Case 1: Both are leaves
    if (A->isLeaf() && B->isLeaf()) {
        // Both zero
        if (A->rank == 0 && B->rank == 0) {
            result->rank = 0;
            result->U = zeroMatrix(A->rows, 0);
            result->V = zeroMatrix(0, A->cols);
            return result;
        }
        
        // Concatenate U matrices and V matrices, then recompress
        Matrix U_combined;
        Matrix V_combined;
        
        if (A->rank > 0) {
            U_combined = A->U;
            V_combined = A->V;
        }
        
        if (B->rank > 0) {
            // Extend U_combined with B->U
            if (U_combined.empty()) {
                U_combined = B->U;
                V_combined = B->V;
            } else {
                for (int i = 0; i < A->rows; ++i) {
                    for (int j = 0; j < B->rank; ++j) {
                        U_combined[i].push_back(B->U[i][j]);
                    }
                }
                for (int i = 0; i < B->rank; ++i) {
                    V_combined.push_back(B->V[i]);
                }
            }
        }
        
        // Recompress: compute full matrix and do SVD
        Matrix dense = matrixMultiply(U_combined, V_combined);
        auto [U_new, S_new, V_new] = svd_decomposition(dense, maxRank, epsilon);
        
        result->rank = static_cast<int>(S_new.size());
        result->U = U_new;
        result->V = V_new;
        
        return result;
    }
    
    // Case 2: Both are internal nodes
    if (!A->isLeaf() && !B->isLeaf()) {
        result->sons.resize(4);
        for (int i = 0; i < 4; ++i) {
            result->sons[i] = hMatrixAdd(A->sons[i], B->sons[i], maxRank, epsilon);
        }
        return result;
    }
    
    // Case 3: Mixed (one leaf, one internal)
    std::shared_ptr<HNode> leaf = A->isLeaf() ? A : B;
    std::shared_ptr<HNode> internal = A->isLeaf() ? B : A;
    
    // Split leaf into 4 sub-blocks by dividing U and V
    int midRow = internal->sons[0]->rows;
    int midCol = internal->sons[0]->cols;
    
    result->sons.resize(4);
    
    if (leaf->rank == 0) {
        // Leaf is zero matrix
        for (int i = 0; i < 4; ++i) {
            result->sons[i] = internal->sons[i];
        }
    } else {
        // Split U into top and bottom
        Matrix U_top = subMatrix(leaf->U, 0, 0, midRow, leaf->rank);
        Matrix U_bottom = subMatrix(leaf->U, midRow, 0, leaf->rows - midRow, leaf->rank);
        
        // Split V into left and right
        Matrix V_left = subMatrix(leaf->V, 0, 0, leaf->rank, midCol);
        Matrix V_right = subMatrix(leaf->V, 0, midCol, leaf->rank, leaf->cols - midCol);
        
        // Create 4 leaf nodes for the split
        auto createLeafNode = [](const Matrix& U, const Matrix& V, int r, int c) {
            auto node = std::make_shared<HNode>(r, c);
            node->rank = rows(U) > 0 ? cols(U) : 0;
            node->U = U;
            node->V = V;
            return node;
        };
        
        auto leaf_TL = createLeafNode(U_top, V_left, midRow, midCol);
        auto leaf_TR = createLeafNode(U_top, V_right, midRow, leaf->cols - midCol);
        auto leaf_BL = createLeafNode(U_bottom, V_left, leaf->rows - midRow, midCol);
        auto leaf_BR = createLeafNode(U_bottom, V_right, leaf->rows - midRow, leaf->cols - midCol);
        
        result->sons[0] = hMatrixAdd(leaf_TL, internal->sons[0], maxRank, epsilon);
        result->sons[1] = hMatrixAdd(leaf_TR, internal->sons[1], maxRank, epsilon);
        result->sons[2] = hMatrixAdd(leaf_BL, internal->sons[2], maxRank, epsilon);
        result->sons[3] = hMatrixAdd(leaf_BR, internal->sons[3], maxRank, epsilon);
    }
    
    return result;
}

// ============================================================================
// H-Matrix Matrix Multiplication (Slides 23-24)
// ============================================================================

std::shared_ptr<HNode> hMatrixMult(const std::shared_ptr<HNode>& A, 
                                    const std::shared_ptr<HNode>& B, 
                                    int maxRank, double epsilon) {
    if (!A || !B || A->cols != B->rows) {
        throw std::runtime_error("Invalid matrix dimensions for multiplication");
    }
    
    auto result = std::make_shared<HNode>(A->rows, B->cols);
    
    // Case 1: Both are leaves
    if (A->isLeaf() && B->isLeaf()) {
        if (A->rank == 0 || B->rank == 0) {
            result->rank = 0;
            result->U = zeroMatrix(A->rows, 0);
            result->V = zeroMatrix(0, B->cols);
            return result;
        }
        
        // Multiply: (U_A * V_A) * (U_B * V_B) = U_A * (V_A * U_B) * V_B
        // Middle part: V_A * U_B is (rank_A × cols_A) * (rows_B × rank_B)
        // Since cols_A = rows_B, this gives (rank_A × rank_B)
        Matrix middle = matrixMultiply(A->V, B->U);
        
        // Now: U_A * middle gives us new U
        Matrix U_result = matrixMultiply(A->U, middle);
        Matrix V_result = B->V;
        
        // Recompress
        Matrix dense = matrixMultiply(U_result, V_result);
        auto [U_new, S_new, V_new] = svd_decomposition(dense, maxRank, epsilon);
        
        result->rank = static_cast<int>(S_new.size());
        result->U = U_new;
        result->V = V_new;
        
        return result;
    }
    
    // Case 2: Both are internal nodes (block multiplication)
    if (!A->isLeaf() && !B->isLeaf()) {
        // Result has 4 sons:
        // [A1 A2] * [B1 B2] = [A1*B1 + A2*B3,  A1*B2 + A2*B4]
        // [A3 A4]   [B3 B4]   [A3*B1 + A4*B3,  A3*B2 + A4*B4]
        
        result->sons.resize(4);
        
        // Top-left: A1*B1 + A2*B3
        auto temp1 = hMatrixMult(A->sons[0], B->sons[0], maxRank, epsilon);
        auto temp2 = hMatrixMult(A->sons[1], B->sons[2], maxRank, epsilon);
        result->sons[0] = hMatrixAdd(temp1, temp2, maxRank, epsilon);
        
        // Top-right: A1*B2 + A2*B4
        temp1 = hMatrixMult(A->sons[0], B->sons[1], maxRank, epsilon);
        temp2 = hMatrixMult(A->sons[1], B->sons[3], maxRank, epsilon);
        result->sons[1] = hMatrixAdd(temp1, temp2, maxRank, epsilon);
        
        // Bottom-left: A3*B1 + A4*B3
        temp1 = hMatrixMult(A->sons[2], B->sons[0], maxRank, epsilon);
        temp2 = hMatrixMult(A->sons[3], B->sons[2], maxRank, epsilon);
        result->sons[2] = hMatrixAdd(temp1, temp2, maxRank, epsilon);
        
        // Bottom-right: A3*B2 + A4*B4
        temp1 = hMatrixMult(A->sons[2], B->sons[1], maxRank, epsilon);
        temp2 = hMatrixMult(A->sons[3], B->sons[3], maxRank, epsilon);
        result->sons[3] = hMatrixAdd(temp1, temp2, maxRank, epsilon);
        
        return result;
    }
    
    // Case 3: Mixed (one leaf, one internal)
    // Convert leaf to internal structure by splitting, then recurse
    std::shared_ptr<HNode> leaf = A->isLeaf() ? A : B;
    std::shared_ptr<HNode> internal = A->isLeaf() ? B : A;
    bool leafIsA = A->isLeaf();
    
    // Split leaf into 4 blocks
    std::shared_ptr<HNode> leafAsInternal = std::make_shared<HNode>(leaf->rows, leaf->cols);
    leafAsInternal->sons.resize(4);
    
    if (leafIsA) {
        int midRow = internal->sons[0]->rows;
        int midCol = leaf->cols / 2;
        
        Matrix U_top = subMatrix(leaf->U, 0, 0, midRow, leaf->rank);
        Matrix U_bottom = subMatrix(leaf->U, midRow, 0, leaf->rows - midRow, leaf->rank);
        Matrix V_left = subMatrix(leaf->V, 0, 0, leaf->rank, midCol);
        Matrix V_right = subMatrix(leaf->V, 0, midCol, leaf->rank, leaf->cols - midCol);
        
        auto createLeaf = [](const Matrix& U, const Matrix& V, int r, int c) {
            auto node = std::make_shared<HNode>(r, c);
            node->rank = cols(U);
            node->U = U;
            node->V = V;
            return node;
        };
        
        leafAsInternal->sons[0] = createLeaf(U_top, V_left, midRow, midCol);
        leafAsInternal->sons[1] = createLeaf(U_top, V_right, midRow, leaf->cols - midCol);
        leafAsInternal->sons[2] = createLeaf(U_bottom, V_left, leaf->rows - midRow, midCol);
        leafAsInternal->sons[3] = createLeaf(U_bottom, V_right, leaf->rows - midRow, leaf->cols - midCol);
        
        return hMatrixMult(leafAsInternal, internal, maxRank, epsilon);
    } else {
        int midRow = leaf->rows / 2;
        int midCol = internal->sons[0]->cols;
        
        Matrix U_top = subMatrix(leaf->U, 0, 0, midRow, leaf->rank);
        Matrix U_bottom = subMatrix(leaf->U, midRow, 0, leaf->rows - midRow, leaf->rank);
        Matrix V_left = subMatrix(leaf->V, 0, 0, leaf->rank, midCol);
        Matrix V_right = subMatrix(leaf->V, 0, midCol, leaf->rank, leaf->cols - midCol);
        
        auto createLeaf = [](const Matrix& U, const Matrix& V, int r, int c) {
            auto node = std::make_shared<HNode>(r, c);
            node->rank = cols(U);
            node->U = U;
            node->V = V;
            return node;
        };
        
        leafAsInternal->sons[0] = createLeaf(U_top, V_left, midRow, midCol);
        leafAsInternal->sons[1] = createLeaf(U_top, V_right, midRow, leaf->cols - midCol);
        leafAsInternal->sons[2] = createLeaf(U_bottom, V_left, leaf->rows - midRow, midCol);
        leafAsInternal->sons[3] = createLeaf(U_bottom, V_right, leaf->rows - midRow, leaf->cols - midCol);
        
        return hMatrixMult(internal, leafAsInternal, maxRank, epsilon);
    }
}

// ============================================================================
// Decompression and Utilities
// ============================================================================

Matrix hMatrixToDense(const std::shared_ptr<HNode>& H) {
    if (!H || H->rows == 0 || H->cols == 0) {
        return zeroMatrix(H ? H->rows : 0, H ? H->cols : 0);
    }
    
    if (H->isLeaf()) {
        if (H->rank == 0) {
            return zeroMatrix(H->rows, H->cols);
        }
        return matrixMultiply(H->U, H->V);
    }
    
    // Recursively decompress children
    Matrix A11 = hMatrixToDense(H->sons[0]);
    Matrix A12 = hMatrixToDense(H->sons[1]);
    Matrix A21 = hMatrixToDense(H->sons[2]);
    Matrix A22 = hMatrixToDense(H->sons[3]);
    
    return combineBlocks(A11, A12, A21, A22);
}

void drawHMatrix(const std::shared_ptr<HNode>& H, Matrix& visualization, 
                 int rowOffset, int colOffset) {
    if (!H) return;
    
    if (H->isLeaf()) {
        // Mark leaf region with value based on rank
        double value = (H->rank > 0) ? 0.0 : 1.0;
        for (int i = 0; i < H->rows; ++i) {
            for (int j = 0; j < H->cols; ++j) {
                visualization[rowOffset + i][colOffset + j] = value;
            }
        }
        
        // Draw rank indicator
        int markSize = std::min(H->rank, std::min(H->rows, H->cols));
        for (int i = 0; i < markSize; ++i) {
            if (rowOffset + i < rows(visualization) && colOffset + i < cols(visualization)) {
                visualization[rowOffset + i][colOffset + i] = 0.5;
            }
        }
    } else {
        // Recurse on children
        int row1 = rowOffset;
        int row2 = rowOffset + H->sons[0]->rows;
        int col1 = colOffset;
        int col2 = colOffset + H->sons[0]->cols;
        
        drawHMatrix(H->sons[0], visualization, row1, col1);
        drawHMatrix(H->sons[1], visualization, row1, col2);
        drawHMatrix(H->sons[2], visualization, row2, col1);
        drawHMatrix(H->sons[3], visualization, row2, col2);
    }
}

Matrix createVisualization(const std::shared_ptr<HNode>& H) {
    if (!H) return Matrix();
    
    Matrix vis = zeroMatrix(H->rows, H->cols);
    for (int i = 0; i < H->rows; ++i) {
        for (int j = 0; j < H->cols; ++j) {
            vis[i][j] = 1.0;
        }
    }
    
    drawHMatrix(H, vis, 0, 0);
    return vis;
}
