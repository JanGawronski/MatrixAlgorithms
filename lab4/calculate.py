"""
Kalkulator mnożenia macierzy przez wektor [1,2,3,4,5,6,7,8]
"""

import numpy as np

# Przykładowa macierz 8x8 (struktura z tablicy)
A = np.array([
    [4, 1, 0, 0, 1, 0, 0, 0],
    [1, 4, 1, 0, 0, 1, 0, 0],
    [0, 1, 4, 1, 0, 0, 1, 0],
    [0, 0, 1, 4, 0, 0, 0, 1],
    [1, 0, 0, 0, 4, 1, 0, 0],
    [0, 1, 0, 0, 1, 4, 1, 0],
    [0, 0, 1, 0, 0, 1, 4, 1],
    [0, 0, 0, 1, 0, 0, 1, 4]
])

# Wektor [1,2,3,4,5,6,7,8]
x = np.array([1, 2, 3, 4, 5, 6, 7, 8])

print("=== Mnożenie macierzy przez wektor [1,2,3,4,5,6,7,8] ===\n")

print("Macierz A:")
print(A)
print()

print("Wektor x:")
print(x)
print()

# Mnożenie
result = A @ x

print("Wynik A * x:")
print(result)
print()

print("Rozpisane:")
for i, val in enumerate(result):
    print(f"y[{i}] = {val}")
print()

print("Wektor wynikowy:")
print(f"[{', '.join(map(str, result))}]")
print()

# Sprawdzenie ręczne dla pierwszego elementu
print("Sprawdzenie pierwszego elementu (y[0]):")
print(f"y[0] = 4*1 + 1*2 + 0*3 + 0*4 + 1*5 + 0*6 + 0*7 + 0*8")
print(f"     = 4 + 2 + 0 + 0 + 5 + 0 + 0 + 0")
print(f"     = {4*1 + 1*2 + 0*3 + 0*4 + 1*5 + 0*6 + 0*7 + 0*8}")
