#include "SupportFunctions.h"
#include "Compression.h"
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image.h"
#include "stb_image_write.h"
#include <iostream>

int main(int argc, char** argv) {
    int rank = 4;
    float epsilon = 1.0;

    if (argc >= 2)
        rank = std::stoi(argv[1]);
    
    if (argc >= 3)
        epsilon = std::stof(argv[2]);

    auto [R, G, B] = loadImageRGB("doge.png");

    auto [U_R, D_R, V_R] = svd_decomposition(R, size(R));
    auto [U_G, D_G, V_G] = svd_decomposition(G, size(G));
    auto [U_B, D_B, V_B] = svd_decomposition(B, size(B));


    for (double val : D_R) std::cout << val << " ";
    std::cout << "\n";
    for (double val : D_G) std::cout << val << " ";
    std::cout << "\n";
    for (double val : D_B) std::cout << val << " ";
    std::cout << "\n";

    TreeNode* rTree = createTree(R, rank, D_R[(D_R.size() - 1) * epsilon]);
    TreeNode* gTree = createTree(G, rank, D_G[(D_G.size() - 1) * epsilon]);
    TreeNode* bTree = createTree(B, rank, D_B[(D_B.size() - 1) * epsilon]);

    Matrix R_comp = drawCompression(rTree, rows(R), cols(R));
    Matrix G_comp = drawCompression(gTree, rows(G), cols(G));
    Matrix B_comp = drawCompression(bTree, rows(B), cols(B));

    Matrix R_rec = reconstructFromTree(rTree);
    Matrix G_rec = reconstructFromTree(gTree);
    Matrix B_rec = reconstructFromTree(bTree);

    saveImageRGB("output/doge_reconstructed_r=" + std::to_string(rank) + "_e=" + std::to_string(epsilon) + ".png", R_rec, G_rec, B_rec);

    saveImageRGB("output/doge_R_reconstructed_r=" + std::to_string(rank) + "_e=" + std::to_string(epsilon) + ".png", R_rec, zeroMatrix(rows(R), cols(R)), zeroMatrix(rows(R), cols(R)));
    saveImageRGB("output/doge_G_reconstructed_r=" + std::to_string(rank) + "_e=" + std::to_string(epsilon) + ".png", zeroMatrix(rows(G), cols(G)), G_rec, zeroMatrix(rows(G), cols(G)));
    saveImageRGB("output/doge_B_reconstructed_r=" + std::to_string(rank) + "_e=" + std::to_string(epsilon) + ".png", zeroMatrix(rows(B), cols(B)), zeroMatrix(rows(B), cols(B)), B_rec);

    saveImageGray("output/doge_R_compressed_r=" + std::to_string(rank) + "_e=" + std::to_string(epsilon) + ".png", R_comp);
    saveImageGray("output/doge_G_compressed_r=" + std::to_string(rank) + "_e=" + std::to_string(epsilon) + ".png", G_comp);
    saveImageGray("output/doge_B_compressed_r=" + std::to_string(rank) + "_e=" + std::to_string(epsilon) + ".png", B_comp);
    return 0;
}
