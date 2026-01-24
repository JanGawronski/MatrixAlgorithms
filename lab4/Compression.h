#pragma once
#include "HMatrix.h"

// Tree node structure for hierarchical compression (from lab3)
// Compatible with HNode structure from HMatrix.h
struct TreeNode {
    Vector singularValues;
    Matrix U;
    Matrix V;
    TreeNode* topLeft = nullptr;
    TreeNode* topRight = nullptr;
    TreeNode* bottomLeft = nullptr;
    TreeNode* bottomRight = nullptr;
};

// Compression tree functions (from lab3)
TreeNode* createTree(const Matrix& A, int rank, double epsilon);
Matrix reconstructFromTree(TreeNode* node);

// Visualization functions (from lab3)
void drawCompression(TreeNode* node, Matrix& matrix, int x0, int x1, int y0, int y1);
Matrix drawCompression(TreeNode* node, int width, int height);

// Conversion functions between TreeNode and HNode
std::shared_ptr<HNode> treeNodeToHNode(TreeNode* tree);
TreeNode* hNodeToTreeNode(const std::shared_ptr<HNode>& h);

// Visualization export functions (saves as PNG images)
bool saveHMatrixVisualizationPNG(const std::shared_ptr<HNode>& H, const std::string& filename);
bool saveTreeVisualizationPNG(TreeNode* tree, int width, int height, const std::string& filename);
