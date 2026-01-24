#include "HMatrix.h"
#include <iostream>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <random>

using namespace std::chrono;

// Save matrix to file for visualization
void saveMatrixToFile(const Matrix& M, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }
    
    for (const auto& row : M) {
        for (size_t i = 0; i < row.size(); ++i) {
            file << row[i];
            if (i < row.size() - 1) file << " ";
        }
        file << "\n";
    }
    file.close();
}

// Save timing results
void saveTimingResults(const std::string& filename, 
                       const std::vector<int>& sizes,
                       const std::vector<double>& times) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }
    
    file << "# Size Time(ms)\n";
    for (size_t i = 0; i < sizes.size(); ++i) {
        file << sizes[i] << " " << times[i] << "\n";
    }
    file.close();
}

int main() {
    std::cout << "=== H-Matrix Implementation - Lab 4 ===" << std::endl;
    std::cout << std::endl;
    
    // Parameters
    std::vector<int> k_values = {2, 3, 4};
    int maxRank = 8;
    double epsilon = 1e-6;
    
    std::vector<int> sizes;
    std::vector<double> vecMultTimes;
    std::vector<double> matMultTimes;
    std::vector<double> vecErrors;
    std::vector<double> matErrors;
    
    for (int k : k_values) {
        int n = 1 << (3 * k);  // 2^(3k)
        sizes.push_back(n);
        
        std::cout << "\n========================================" << std::endl;
        std::cout << "k = " << k << ", Matrix size = " << n << " x " << n << std::endl;
        std::cout << "========================================" << std::endl;
        
        // Step 1: Generate 3D grid matrix
        std::cout << "\n[1/7] Generating 3D grid matrix..." << std::flush;
        auto start = high_resolution_clock::now();
        Matrix A = generate3DGridMatrix(k);
        auto end = high_resolution_clock::now();
        std::cout << " done (" << duration_cast<milliseconds>(end - start).count() << " ms)" << std::endl;
        
        // Step 2: Build H-Matrix
        std::cout << "[2/7] Building H-Matrix (rank=" << maxRank << ", epsilon=" << epsilon << ")..." << std::flush;
        start = high_resolution_clock::now();
        auto H = buildHMatrix(A, maxRank, epsilon);
        end = high_resolution_clock::now();
        std::cout << " done (" << duration_cast<milliseconds>(end - start).count() << " ms)" << std::endl;
        
        // Step 3: Visualize H-Matrix structure
        std::cout << "[3/7] Creating H-Matrix visualization..." << std::flush;
        Matrix vis = createVisualization(H);
        std::string visFilename = "hmatrix_structure_k" + std::to_string(k) + ".txt";
        saveMatrixToFile(vis, visFilename);
        std::cout << " saved to " << visFilename << std::endl;
        
        // Step 4: Matrix-Vector multiplication
        std::cout << "[4/7] Testing H-Matrix * vector..." << std::flush;
        
        // Generate random vector
        Vector x(n);
        std::mt19937 gen(42);
        std::uniform_real_distribution<> dis(0.0, 1.0);
        for (int i = 0; i < n; ++i) x[i] = dis(gen);
        
        // Time H-Matrix multiplication
        start = high_resolution_clock::now();
        Vector Hx = hMatrixVectorMult(H, x);
        end = high_resolution_clock::now();
        double hvTime = duration_cast<microseconds>(end - start).count() / 1000.0;
        vecMultTimes.push_back(hvTime);
        std::cout << " done (" << hvTime << " ms)" << std::endl;
        
        // Compute dense multiplication for comparison
        std::cout << "    Computing dense matrix * vector for comparison..." << std::flush;
        start = high_resolution_clock::now();
        Vector Ax = matrixVectorMult(A, x);
        end = high_resolution_clock::now();
        std::cout << " done (" << duration_cast<milliseconds>(end - start).count() << " ms)" << std::endl;
        
        // Compute error
        double vecError = std::sqrt(vectorError(Ax, Hx));
        vecErrors.push_back(vecError);
        std::cout << "    Error ||Ax - Hx||_2 = " << std::scientific << vecError << std::endl;
        
        // Step 5: Matrix-Matrix multiplication (A^2)
        std::cout << "[5/7] Testing H-Matrix * H-Matrix (squaring)..." << std::flush;
        
        start = high_resolution_clock::now();
        auto H2 = hMatrixMult(H, H, maxRank, epsilon);
        end = high_resolution_clock::now();
        double mmTime = duration_cast<milliseconds>(end - start).count();
        matMultTimes.push_back(mmTime);
        std::cout << " done (" << mmTime << " ms)" << std::endl;
        
        // Step 6: Compute dense A^2 for comparison (only for smaller sizes)
        if (n <= 512) {
            std::cout << "[6/7] Computing dense A^2 for comparison..." << std::flush;
            start = high_resolution_clock::now();
            Matrix A2_dense = matrixMultiply(A, A);
            end = high_resolution_clock::now();
            std::cout << " done (" << duration_cast<milliseconds>(end - start).count() << " ms)" << std::endl;
            
            std::cout << "[7/7] Decompressing H^2 to dense form..." << std::flush;
            start = high_resolution_clock::now();
            Matrix H2_dense = hMatrixToDense(H2);
            end = high_resolution_clock::now();
            std::cout << " done (" << duration_cast<milliseconds>(end - start).count() << " ms)" << std::endl;
            
            // Compute Frobenius norm error
            double sumSq = 0.0;
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    double diff = A2_dense[i][j] - H2_dense[i][j];
                    sumSq += diff * diff;
                }
            }
            double matError = std::sqrt(sumSq);
            matErrors.push_back(matError);
            std::cout << "    Error ||A^2 - H^2||_F = " << std::scientific << matError << std::endl;
        } else {
            std::cout << "[6/7] Skipping dense computation (matrix too large)" << std::endl;
            std::cout << "[7/7] Skipping decompression (matrix too large)" << std::endl;
            matErrors.push_back(-1.0);  // Mark as not computed
        }
        
        std::cout << std::endl;
    }
    
    // Save timing results
    std::cout << "\n=== Saving Results ===" << std::endl;
    saveTimingResults("timing_vector_mult.txt", sizes, vecMultTimes);
    std::cout << "Vector multiplication times saved to timing_vector_mult.txt" << std::endl;
    
    saveTimingResults("timing_matrix_mult.txt", sizes, matMultTimes);
    std::cout << "Matrix multiplication times saved to timing_matrix_mult.txt" << std::endl;
    
    // Save errors
    std::ofstream errFile("errors.txt");
    errFile << "# Size VectorError MatrixError\n";
    for (size_t i = 0; i < sizes.size(); ++i) {
        errFile << sizes[i] << " " << std::scientific << vecErrors[i] << " " << matErrors[i] << "\n";
    }
    errFile.close();
    std::cout << "Errors saved to errors.txt" << std::endl;
    
    // Summary table
    std::cout << "\n=== Summary ===" << std::endl;
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "┌───────┬──────────┬───────────────┬───────────────┬─────────────┬─────────────┐" << std::endl;
    std::cout << "│   k   │   Size   │  H*v time(ms) │  H*H time(ms) │ Vector Error│ Matrix Error│" << std::endl;
    std::cout << "├───────┼──────────┼───────────────┼───────────────┼─────────────┼─────────────┤" << std::endl;
    for (size_t i = 0; i < sizes.size(); ++i) {
        std::cout << "│  " << std::setw(4) << k_values[i] << " │ " 
                  << std::setw(8) << sizes[i] << " │ "
                  << std::setw(13) << vecMultTimes[i] << " │ "
                  << std::setw(13) << matMultTimes[i] << " │ "
                  << std::scientific << std::setprecision(2) << std::setw(11) << vecErrors[i] << " │ ";
        if (matErrors[i] >= 0) {
            std::cout << std::setw(11) << matErrors[i];
        } else {
            std::cout << "     N/A    ";
        }
        std::cout << " │" << std::endl;
    }
    std::cout << "└───────┴──────────┴───────────────┴───────────────┴─────────────┴─────────────┘" << std::endl;
    
    // Complexity analysis
    std::cout << "\n=== Complexity Analysis ===" << std::endl;
    std::cout << "For matrix-vector multiplication:" << std::endl;
    if (vecMultTimes.size() >= 2) {
        double ratio1 = vecMultTimes[1] / vecMultTimes[0];
        double size_ratio1 = static_cast<double>(sizes[1]) / sizes[0];
        double beta1 = std::log(ratio1) / std::log(size_ratio1);
        std::cout << "  Time ratio (k=3/k=2): " << ratio1 << std::endl;
        std::cout << "  Size ratio: " << size_ratio1 << std::endl;
        std::cout << "  Estimated β (T = αN^β): " << std::fixed << std::setprecision(3) << beta1 << std::endl;
    }
    if (vecMultTimes.size() >= 3) {
        double ratio2 = vecMultTimes[2] / vecMultTimes[1];
        double size_ratio2 = static_cast<double>(sizes[2]) / sizes[1];
        double beta2 = std::log(ratio2) / std::log(size_ratio2);
        std::cout << "  Time ratio (k=4/k=3): " << ratio2 << std::endl;
        std::cout << "  Size ratio: " << size_ratio2 << std::endl;
        std::cout << "  Estimated β (T = αN^β): " << std::fixed << std::setprecision(3) << beta2 << std::endl;
    }
    
    std::cout << "\nFor matrix-matrix multiplication:" << std::endl;
    if (matMultTimes.size() >= 2) {
        double ratio1 = matMultTimes[1] / matMultTimes[0];
        double size_ratio1 = static_cast<double>(sizes[1]) / sizes[0];
        double beta1 = std::log(ratio1) / std::log(size_ratio1);
        std::cout << "  Time ratio (k=3/k=2): " << ratio1 << std::endl;
        std::cout << "  Size ratio: " << size_ratio1 << std::endl;
        std::cout << "  Estimated β (T = αN^β): " << std::fixed << std::setprecision(3) << beta1 << std::endl;
    }
    if (matMultTimes.size() >= 3) {
        double ratio2 = matMultTimes[2] / matMultTimes[1];
        double size_ratio2 = static_cast<double>(sizes[2]) / sizes[1];
        double beta2 = std::log(ratio2) / std::log(size_ratio2);
        std::cout << "  Time ratio (k=4/k=3): " << ratio2 << std::endl;
        std::cout << "  Size ratio: " << size_ratio2 << std::endl;
        std::cout << "  Estimated β (T = αN^β): " << std::fixed << std::setprecision(3) << beta2 << std::endl;
    }
    
    std::cout << "\n=== Program completed successfully ===" << std::endl;
    
    return 0;
}
