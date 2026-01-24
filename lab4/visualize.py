#!/usr/bin/env python3
"""
Skrypt do wizualizacji wyników H-Matrix
Generuje wykresy z danych timing i errors
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def load_timing_data(filename):
    """Wczytaj dane czasowe z pliku"""
    data = np.loadtxt(filename, comments='#')
    if data.ndim == 1:
        data = data.reshape(1, -1)
    sizes = data[:, 0].astype(int)
    times = data[:, 1]
    return sizes, times

def load_error_data(filename):
    """Wczytaj dane o błędach z pliku"""
    data = np.loadtxt(filename, comments='#')
    if data.ndim == 1:
        data = data.reshape(1, -1)
    sizes = data[:, 0].astype(int)
    vec_errors = data[:, 1]
    mat_errors = data[:, 2]
    return sizes, vec_errors, mat_errors

def fit_power_law(x, y):
    """Dopasuj T = α * N^β"""
    # log(T) = log(α) + β * log(N)
    log_x = np.log(x)
    log_y = np.log(y)
    
    # Least squares fit
    A = np.vstack([log_x, np.ones(len(log_x))]).T
    beta, log_alpha = np.linalg.lstsq(A, log_y, rcond=None)[0]
    alpha = np.exp(log_alpha)
    
    return alpha, beta

def plot_timing_results():
    """Generuj wykresy czasów"""
    
    # Wczytaj dane
    try:
        sizes_vec, times_vec = load_timing_data('timing_vector_mult.txt')
        sizes_mat, times_mat = load_timing_data('timing_matrix_mult.txt')
    except:
        print("Błąd wczytywania danych. Upewnij się, że pliki istnieją.")
        return
    
    # Dopasuj funkcje potęgowe
    if len(sizes_vec) >= 2:
        alpha_vec, beta_vec = fit_power_law(sizes_vec, times_vec)
        fit_vec = alpha_vec * sizes_vec**beta_vec
    
    if len(sizes_mat) >= 2:
        alpha_mat, beta_mat = fit_power_law(sizes_mat, times_mat)
        fit_mat = alpha_mat * sizes_mat**beta_mat
    
    # Wykres dla mnożenia macierz-wektor
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    ax1.plot(sizes_vec, times_vec, 'bo-', label='Dane eksperymentalne', markersize=8, linewidth=2)
    if len(sizes_vec) >= 2:
        ax1.plot(sizes_vec, fit_vec, 'r--', 
                label=f'Dopasowanie: T = {alpha_vec:.2e} × N^{beta_vec:.2f}', linewidth=2)
    ax1.set_xlabel('Rozmiar macierzy N', fontsize=12)
    ax1.set_ylabel('Czas [ms]', fontsize=12)
    ax1.set_title('Mnożenie H-Matrix × wektor', fontsize=14, fontweight='bold')
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)
    
    # Wykres dla mnożenia macierz-macierz
    ax2.plot(sizes_mat, times_mat, 'rs-', label='Dane eksperymentalne', markersize=8, linewidth=2)
    if len(sizes_mat) >= 2:
        ax2.plot(sizes_mat, fit_mat, 'b--', 
                label=f'Dopasowanie: T = {alpha_mat:.2e} × N^{beta_mat:.2f}', linewidth=2)
    ax2.set_xlabel('Rozmiar macierzy N', fontsize=12)
    ax2.set_ylabel('Czas [ms]', fontsize=12)
    ax2.set_title('Mnożenie H-Matrix × H-Matrix', fontsize=14, fontweight='bold')
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3)
    ax2.set_yscale('log')
    
    plt.tight_layout()
    plt.savefig('timing_plots.png', dpi=300, bbox_inches='tight')
    print(f"Zapisano wykres do: timing_plots.png")
    
    # Wydrukuj wyniki
    print("\n=== Analiza złożoności ===")
    if len(sizes_vec) >= 2:
        print(f"\nMnożenie macierz-wektor:")
        print(f"  T = {alpha_vec:.2e} × N^{beta_vec:.3f}")
        print(f"  β = {beta_vec:.3f} (oczekiwane: 1.0-1.5 dla O(N log N))")
    
    if len(sizes_mat) >= 2:
        print(f"\nMnożenie macierz-macierz:")
        print(f"  T = {alpha_mat:.2e} × N^{beta_mat:.3f}")
        print(f"  β = {beta_mat:.3f} (oczekiwane: 2.0-2.5 dla O(N²))")
    
    plt.show()

def plot_errors():
    """Generuj wykresy błędów"""
    
    try:
        sizes, vec_errors, mat_errors = load_error_data('errors.txt')
    except:
        print("Błąd wczytywania danych o błędach.")
        return
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    # Błąd macierz-wektor
    ax1.semilogy(sizes, vec_errors, 'go-', markersize=8, linewidth=2)
    ax1.set_xlabel('Rozmiar macierzy N', fontsize=12)
    ax1.set_ylabel('||Ax - Hx||₂', fontsize=12)
    ax1.set_title('Błąd mnożenia macierz-wektor', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    
    # Błąd macierz-macierz (pomijając -1)
    valid_mat = mat_errors > 0
    if np.any(valid_mat):
        ax2.semilogy(sizes[valid_mat], mat_errors[valid_mat], 'mo-', markersize=8, linewidth=2)
        ax2.set_xlabel('Rozmiar macierzy N', fontsize=12)
        ax2.set_ylabel('||A² - H²||_F', fontsize=12)
        ax2.set_title('Błąd mnożenia macierz-macierz', fontsize=14, fontweight='bold')
        ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('error_plots.png', dpi=300, bbox_inches='tight')
    print(f"Zapisano wykres błędów do: error_plots.png")
    plt.show()

def visualize_hmatrix_structure(k):
    """Wizualizuj strukturę H-macierzy"""
    filename = f'hmatrix_structure_k{k}.txt'
    
    try:
        data = np.loadtxt(filename)
    except:
        print(f"Nie znaleziono pliku: {filename}")
        return
    
    plt.figure(figsize=(8, 8))
    plt.imshow(data, cmap='RdYlGn', interpolation='nearest')
    plt.title(f'Struktura H-macierzy dla k={k} (N={data.shape[0]})', 
              fontsize=14, fontweight='bold')
    plt.colorbar(label='0=skompresowane, 0.5=rank marker, 1=podzielone')
    plt.xlabel('Kolumna', fontsize=12)
    plt.ylabel('Wiersz', fontsize=12)
    
    filename_out = f'hmatrix_structure_k{k}.png'
    plt.savefig(filename_out, dpi=300, bbox_inches='tight')
    print(f"Zapisano wizualizację do: {filename_out}")
    plt.show()

def main():
    print("=== Wizualizacja wyników H-Matrix ===\n")
    
    # 1. Wykresy czasowe
    print("[1/3] Generowanie wykresów czasowych...")
    plot_timing_results()
    
    # 2. Wykresy błędów
    print("\n[2/3] Generowanie wykresów błędów...")
    plot_errors()
    
    # 3. Wizualizacje struktur
    print("\n[3/3] Generowanie wizualizacji struktur...")
    for k in [2, 3, 4]:
        visualize_hmatrix_structure(k)
    
    print("\n=== Zakończono ===")

if __name__ == '__main__':
    main()
