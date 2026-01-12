#pragma once
#include "SupportFunctions.h"

struct TreeNode {
    Vector singularValues;
    Matrix U;
    Matrix V;
    TreeNode* topLeft = nullptr;
    TreeNode* topRight = nullptr;
    TreeNode* bottomLeft = nullptr;
    TreeNode* bottomRight = nullptr;
};


TreeNode* createTree(const Matrix &A, int rank, double epsilon);
Matrix reconstructFromTree(TreeNode* node);
void drawCompression(TreeNode* node, Matrix &matrix, int x0, int x1, int y0, int y1);
Matrix drawCompression(TreeNode* node, int width, int height);