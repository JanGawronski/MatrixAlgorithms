import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

# Read Binet data
binet_data = {'N': [], 'time': [], 'adds': [], 'subs': [], 'muls': [], 'divs': [], 'memory': []}
with open('wynikBinet.txt', 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.strip().split()
        if len(parts) >= 8:
            binet_data['N'].append(int(parts[0]))
            binet_data['time'].append(float(parts[1]))
            binet_data['adds'].append(int(parts[2]))
            binet_data['subs'].append(int(parts[3]))
            binet_data['muls'].append(int(parts[4]))
            binet_data['divs'].append(int(parts[5]))
            binet_data['memory'].append(int(parts[6]))

# Read Strassen data
strassen_data = {'N': [], 'time': [], 'adds': [], 'subs': [], 'muls': [], 'divs': [], 'memory': []}
with open('wynikStrassen.txt', 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.strip().split()
        if len(parts) >= 8:
            strassen_data['N'].append(int(parts[0]))
            strassen_data['time'].append(float(parts[1]))
            strassen_data['adds'].append(int(parts[2]))
            strassen_data['subs'].append(int(parts[3]))
            strassen_data['muls'].append(int(parts[4]))
            strassen_data['divs'].append(int(parts[5]))
            strassen_data['memory'].append(int(parts[6]))

# Calculate total floating point operations
binet_ops = [binet_data['adds'][i] + binet_data['subs'][i] + binet_data['muls'][i] + binet_data['divs'][i] 
             for i in range(len(binet_data['N']))]
strassen_ops = [strassen_data['adds'][i] + strassen_data['subs'][i] + strassen_data['muls'][i] + strassen_data['divs'][i] 
                for i in range(len(strassen_data['N']))]

# Set up the style
plt.style.use('seaborn-v0_8-whitegrid')
fig_size = (10, 6)

# 1. Wykres czasu działania
plt.figure(figsize=fig_size)
plt.plot(binet_data['N'], binet_data['time'], 'o-', label='Binet (rekurencyjny)', linewidth=2, markersize=5, color='#2E86AB')
plt.plot(strassen_data['N'], strassen_data['time'], 's-', label='Strassen', linewidth=2, markersize=5, color='#A23B72')
plt.xlabel('Rozmiar macierzy N', fontsize=13, fontweight='bold')
plt.ylabel('Czas wykonania [s]', fontsize=13, fontweight='bold')
plt.title('Porównanie czasu wykonania algorytmów mnożenia macierzy', fontsize=15, fontweight='bold', pad=15)
plt.legend(fontsize=12, loc='upper left')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('wykres_czas.png', dpi=300, bbox_inches='tight')
print("✓ Wygenerowano: wykres_czas.png")
plt.close()

# 2. Wykres czasu działania (skala logarytmiczna)
plt.figure(figsize=fig_size)
plt.semilogy(binet_data['N'], binet_data['time'], 'o-', label='Binet (rekurencyjny)', linewidth=2, markersize=5, color='#2E86AB')
plt.semilogy(strassen_data['N'], strassen_data['time'], 's-', label='Strassen', linewidth=2, markersize=5, color='#A23B72')
plt.xlabel('Rozmiar macierzy N', fontsize=13, fontweight='bold')
plt.ylabel('Czas wykonania [s] (skala logarytmiczna)', fontsize=13, fontweight='bold')
plt.title('Czas wykonania - skala logarytmiczna', fontsize=15, fontweight='bold', pad=15)
plt.legend(fontsize=12, loc='upper left')
plt.grid(True, alpha=0.3, which='both')
plt.tight_layout()
plt.savefig('wykres_czas_log.png', dpi=300, bbox_inches='tight')
print("✓ Wygenerowano: wykres_czas_log.png")
plt.close()

# 3. Wykres liczby operacji zmiennoprzecinkowych
plt.figure(figsize=fig_size)
plt.plot(binet_data['N'], binet_ops, 'o-', label='Binet (rekurencyjny)', linewidth=2, markersize=5, color='#2E86AB')
plt.plot(strassen_data['N'], strassen_ops, 's-', label='Strassen', linewidth=2, markersize=5, color='#A23B72')
plt.xlabel('Rozmiar macierzy N', fontsize=13, fontweight='bold')
plt.ylabel('Liczba operacji zmiennoprzecinkowych', fontsize=13, fontweight='bold')
plt.title('Liczba operacji zmiennoprzecinkowych', fontsize=15, fontweight='bold', pad=15)
plt.legend(fontsize=12, loc='upper left')
plt.grid(True, alpha=0.3)
plt.ticklabel_format(style='plain', axis='y')
plt.tight_layout()
plt.savefig('wykres_operacje.png', dpi=300, bbox_inches='tight')
print("✓ Wygenerowano: wykres_operacje.png")
plt.close()

# 4. Wykres operacji (skala logarytmiczna)
plt.figure(figsize=fig_size)
plt.loglog(binet_data['N'], binet_ops, 'o-', label='Binet (rekurencyjny)', linewidth=2, markersize=5, color='#2E86AB')
plt.loglog(strassen_data['N'], strassen_ops, 's-', label='Strassen', linewidth=2, markersize=5, color='#A23B72')
plt.xlabel('Rozmiar macierzy N (skala logarytmiczna)', fontsize=13, fontweight='bold')
plt.ylabel('Liczba operacji (skala logarytmiczna)', fontsize=13, fontweight='bold')
plt.title('Liczba operacji - skala logarytmiczna', fontsize=15, fontweight='bold', pad=15)
plt.legend(fontsize=12, loc='upper left')
plt.grid(True, alpha=0.3, which='both')
plt.tight_layout()
plt.savefig('wykres_operacje_log.png', dpi=300, bbox_inches='tight')
print("✓ Wygenerowano: wykres_operacje_log.png")
plt.close()

# 5. Wykres zużycia pamięci
plt.figure(figsize=fig_size)
plt.plot(binet_data['N'], [m/1024 for m in binet_data['memory']], 'o-', label='Binet (rekurencyjny)', linewidth=2, markersize=5, color='#2E86AB')
plt.plot(strassen_data['N'], [m/1024 for m in strassen_data['memory']], 's-', label='Strassen', linewidth=2, markersize=5, color='#A23B72')
plt.xlabel('Rozmiar macierzy N', fontsize=13, fontweight='bold')
plt.ylabel('Zużycie pamięci [KB]', fontsize=13, fontweight='bold')
plt.title('Szczytowe zużycie pamięci', fontsize=15, fontweight='bold', pad=15)
plt.legend(fontsize=12, loc='upper left')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('wykres_pamiec.png', dpi=300, bbox_inches='tight')
print("✓ Wygenerowano: wykres_pamiec.png")
plt.close()

# 6. Wykres pamięci (skala logarytmiczna)
plt.figure(figsize=fig_size)
plt.loglog(binet_data['N'], [m/1024 for m in binet_data['memory']], 'o-', label='Binet (rekurencyjny)', linewidth=2, markersize=5, color='#2E86AB')
plt.loglog(strassen_data['N'], [m/1024 for m in strassen_data['memory']], 's-', label='Strassen', linewidth=2, markersize=5, color='#A23B72')
plt.xlabel('Rozmiar macierzy N (skala logarytmiczna)', fontsize=13, fontweight='bold')
plt.ylabel('Zużycie pamięci [KB] (skala logarytmiczna)', fontsize=13, fontweight='bold')
plt.title('Zużycie pamięci - skala logarytmiczna', fontsize=15, fontweight='bold', pad=15)
plt.legend(fontsize=12, loc='upper left')
plt.grid(True, alpha=0.3, which='both')
plt.tight_layout()
plt.savefig('wykres_pamiec_log.png', dpi=300, bbox_inches='tight')
print("✓ Wygenerowano: wykres_pamiec_log.png")
plt.close()

# 7. Analiza złożoności - dopasowanie krzywej
def power_law(x, a, b):
    return a * np.power(x, b)

# Skip first few points for better fit
start_idx = 5
popt_binet, _ = curve_fit(power_law, binet_data['N'][start_idx:], binet_ops[start_idx:], p0=[1, 3])
popt_strassen, _ = curve_fit(power_law, strassen_data['N'][start_idx:], strassen_ops[start_idx:], p0=[1, 2.8])

plt.figure(figsize=fig_size)
plt.loglog(binet_data['N'], binet_ops, 'o', label='Binet (dane pomiarowe)', alpha=0.6, markersize=6, color='#2E86AB')
plt.loglog(strassen_data['N'], strassen_ops, 's', label='Strassen (dane pomiarowe)', alpha=0.6, markersize=6, color='#A23B72')

# Fitted curves
x_fit = np.array(binet_data['N'])
plt.loglog(x_fit, power_law(x_fit, *popt_binet), '--', 
           label=f'Binet dopasowanie: O(n^{popt_binet[1]:.2f})', linewidth=2.5, color='#2E86AB')
plt.loglog(x_fit, power_law(x_fit, *popt_strassen), '--', 
           label=f'Strassen dopasowanie: O(n^{popt_strassen[1]:.2f})', linewidth=2.5, color='#A23B72')

plt.xlabel('Rozmiar macierzy N (skala logarytmiczna)', fontsize=13, fontweight='bold')
plt.ylabel('Liczba operacji (skala logarytmiczna)', fontsize=13, fontweight='bold')
plt.title('Analiza złożoności obliczeniowej', fontsize=15, fontweight='bold', pad=15)
plt.legend(fontsize=10, loc='upper left')
plt.grid(True, alpha=0.3, which='both')
plt.tight_layout()
plt.savefig('wykres_zlozonosc.png', dpi=300, bbox_inches='tight')
print("✓ Wygenerowano: wykres_zlozonosc.png")
plt.close()

# 8. Porównanie typów operacji dla Binet
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

axes[0, 0].plot(binet_data['N'], binet_data['adds'], 'o-', color='#06D6A0', linewidth=2, markersize=4)
axes[0, 0].set_xlabel('Rozmiar macierzy N', fontsize=11)
axes[0, 0].set_ylabel('Liczba operacji', fontsize=11)
axes[0, 0].set_title('Binet - Dodawania', fontsize=12, fontweight='bold')
axes[0, 0].grid(True, alpha=0.3)

axes[0, 1].plot(binet_data['N'], binet_data['subs'], 'o-', color='#FFD166', linewidth=2, markersize=4)
axes[0, 1].set_xlabel('Rozmiar macierzy N', fontsize=11)
axes[0, 1].set_ylabel('Liczba operacji', fontsize=11)
axes[0, 1].set_title('Binet - Odejmowania', fontsize=12, fontweight='bold')
axes[0, 1].grid(True, alpha=0.3)

axes[1, 0].plot(binet_data['N'], binet_data['muls'], 'o-', color='#EF476F', linewidth=2, markersize=4)
axes[1, 0].set_xlabel('Rozmiar macierzy N', fontsize=11)
axes[1, 0].set_ylabel('Liczba operacji', fontsize=11)
axes[1, 0].set_title('Binet - Mnożenia', fontsize=12, fontweight='bold')
axes[1, 0].grid(True, alpha=0.3)

axes[1, 1].plot(binet_data['N'], binet_data['adds'], 'o-', label='Dodawania', linewidth=2, markersize=4)
axes[1, 1].plot(binet_data['N'], binet_data['subs'], 's-', label='Odejmowania', linewidth=2, markersize=4)
axes[1, 1].plot(binet_data['N'], binet_data['muls'], '^-', label='Mnożenia', linewidth=2, markersize=4)
axes[1, 1].set_xlabel('Rozmiar macierzy N', fontsize=11)
axes[1, 1].set_ylabel('Liczba operacji', fontsize=11)
axes[1, 1].set_title('Binet - Wszystkie typy operacji', fontsize=12, fontweight='bold')
axes[1, 1].legend(fontsize=10)
axes[1, 1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('wykres_operacje_szczegoly_binet.png', dpi=300, bbox_inches='tight')
print("✓ Wygenerowano: wykres_operacje_szczegoly_binet.png")
plt.close()

# 9. Porównanie typów operacji dla Strassen
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

axes[0, 0].plot(strassen_data['N'], strassen_data['adds'], 'o-', color='#06D6A0', linewidth=2, markersize=4)
axes[0, 0].set_xlabel('Rozmiar macierzy N', fontsize=11)
axes[0, 0].set_ylabel('Liczba operacji', fontsize=11)
axes[0, 0].set_title('Strassen - Dodawania', fontsize=12, fontweight='bold')
axes[0, 0].grid(True, alpha=0.3)

axes[0, 1].plot(strassen_data['N'], strassen_data['subs'], 'o-', color='#FFD166', linewidth=2, markersize=4)
axes[0, 1].set_xlabel('Rozmiar macierzy N', fontsize=11)
axes[0, 1].set_ylabel('Liczba operacji', fontsize=11)
axes[0, 1].set_title('Strassen - Odejmowania', fontsize=12, fontweight='bold')
axes[0, 1].grid(True, alpha=0.3)

axes[1, 0].plot(strassen_data['N'], strassen_data['muls'], 'o-', color='#EF476F', linewidth=2, markersize=4)
axes[1, 0].set_xlabel('Rozmiar macierzy N', fontsize=11)
axes[1, 0].set_ylabel('Liczba operacji', fontsize=11)
axes[1, 0].set_title('Strassen - Mnożenia', fontsize=12, fontweight='bold')
axes[1, 0].grid(True, alpha=0.3)

axes[1, 1].plot(strassen_data['N'], strassen_data['adds'], 'o-', label='Dodawania', linewidth=2, markersize=4)
axes[1, 1].plot(strassen_data['N'], strassen_data['subs'], 's-', label='Odejmowania', linewidth=2, markersize=4)
axes[1, 1].plot(strassen_data['N'], strassen_data['muls'], '^-', label='Mnożenia', linewidth=2, markersize=4)
axes[1, 1].set_xlabel('Rozmiar macierzy N', fontsize=11)
axes[1, 1].set_ylabel('Liczba operacji', fontsize=11)
axes[1, 1].set_title('Strassen - Wszystkie typy operacji', fontsize=12, fontweight='bold')
axes[1, 1].legend(fontsize=10)
axes[1, 1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('wykres_operacje_szczegoly_strassen.png', dpi=300, bbox_inches='tight')
print("✓ Wygenerowano: wykres_operacje_szczegoly_strassen.png")
plt.close()

print("\n" + "="*60)
print("PODSUMOWANIE ANALIZY ZŁOŻONOŚCI OBLICZENIOWEJ")
print("="*60)
print(f"\nBinet (rekurencyjny):")
print(f"  Dopasowana złożoność: O(n^{popt_binet[1]:.3f})")
print(f"  Złożoność teoretyczna: O(n^3)")
print(f"  Różnica: {abs(popt_binet[1] - 3.0):.3f}")

print(f"\nStrassen:")
print(f"  Dopasowana złożoność: O(n^{popt_strassen[1]:.3f})")
print(f"  Złożoność teoretyczna: O(n^2.807)")
print(f"  Różnica: {abs(popt_strassen[1] - 2.807):.3f}")

print(f"\nZakres pomiarów:")
print(f"  Rozmiary macierzy: {min(binet_data['N'])} - {max(binet_data['N'])}")
print(f"  Liczba pomiarów: {len(binet_data['N'])}")

print("\n" + "="*60)
print("Wszystkie wykresy zostały wygenerowane pomyślnie!")
print("="*60)