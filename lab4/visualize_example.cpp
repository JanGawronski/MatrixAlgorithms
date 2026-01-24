#include "HMatrix.h"
#include "Compression.h"
#include <iostream>
#include <iomanip>

int main() {
    std::cout << "=== Przykład wizualizacji H-macierzy ===" << std::endl;
    std::cout << std::endl;
    
    // Parametry
    int k = 2;  // Dla macierzy 64x64
    int maxRank = 8;
    double epsilon = 1e-6;
    
    std::cout << "Generowanie macierzy 3D grid dla k=" << k << " (rozmiar: " 
              << (1 << (3*k)) << "x" << (1 << (3*k)) << ")..." << std::endl;
    
    Matrix A = generate3DGridMatrix(k);
    int n = rows(A);
    
    std::cout << "Budowanie H-macierzy (rank=" << maxRank << ", epsilon=" 
              << epsilon << ")..." << std::endl;
    
    auto H = buildHMatrix(A, maxRank, epsilon);
    
    std::cout << "H-macierz zbudowana." << std::endl;
    std::cout << std::endl;
    
    // Metoda 1: Użycie HNode bezpośrednio
    std::cout << "Tworzenie wizualizacji metoda 1 (HNode)..." << std::endl;
    Matrix vis1 = createVisualization(H);
    
    std::string filename1 = "hmatrix_vis_method1.png";
    if (saveImageGray(filename1, vis1)) {
        std::cout << "✓ Zapisano wizualizację do: " << filename1 << std::endl;
    } else {
        std::cout << "✗ Błąd zapisu do: " << filename1 << std::endl;
    }
    
    // Metoda 2: Konwersja do TreeNode (kompatybilność z Lab3)
    std::cout << "\nTworzenie wizualizacji metoda 2 (TreeNode z Lab3)..." << std::endl;
    TreeNode* tree = hNodeToTreeNode(H);
    Matrix vis2 = drawCompression(tree, n, n);
    
    std::string filename2 = "hmatrix_vis_method2.png";
    if (saveImageGray(filename2, vis2)) {
        std::cout << "✓ Zapisano wizualizację do: " << filename2 << std::endl;
    } else {
        std::cout << "✗ Błąd zapisu do: " << filename2 << std::endl;
    }
    
    // Metoda 3: Bezpośrednie użycie createTree z Lab3
    std::cout << "\nTworzenie wizualizacji metoda 3 (createTree z Lab3)..." << std::endl;
    TreeNode* tree2 = createTree(A, maxRank, epsilon);
    
    std::string filename3 = "hmatrix_vis_method3.png";
    if (saveTreeVisualizationPNG(tree2, n, n, filename3)) {
        std::cout << "✓ Zapisano wizualizację do: " << filename3 << std::endl;
    } else {
        std::cout << "✗ Błąd zapisu do: " << filename3 << std::endl;
    }
    
    // Wyświetl fragment wizualizacji w terminalu
    std::cout << "\nFragment wizualizacji (8x8 lewy górny róg):" << std::endl;
    std::cout << "Legenda: 1.0=podzielone, 0.0=skompresowane, 0.5=marker ranku" << std::endl;
    std::cout << std::endl;
    for (int i = 0; i < std::min(8, rows(vis1)); ++i) {
        for (int j = 0; j < std::min(8, cols(vis1)); ++j) {
            if (vis1[i][j] > 0.75) std::cout << "█ ";
            else if (vis1[i][j] > 0.25) std::cout << "▓ ";
            else std::cout << "░ ";
        }
        std::cout << std::endl;
    }
    
    std::cout << "\n=== Zakończono ===" << std::endl;
    std::cout << "\nZapisano 3 pliki PNG z wizualizacjami:" << std::endl;
    std::cout << "  - " << filename1 << " (bezpośrednio z HNode)" << std::endl;
    std::cout << "  - " << filename2 << " (konwersja HNode→TreeNode)" << std::endl;
    std::cout << "  - " << filename3 << " (createTree z Lab3)" << std::endl;
    
    return 0;
}
