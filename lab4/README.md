# Lab 4 - H-Matrix Implementation

## Opis zadania

Implementacja hierarchicznych macierzy (H-Matrix) dla macierzy reprezentujących topologię 3D siatki sześciennej.

## Struktura projektu

```
lab4/
├── HMatrix.h              - Nagłówek z definicjami struktur i funkcji H-macierzy
├── HMatrix.cpp            - Implementacja operacji na H-macierzach
├── SupportFunctions.cpp   - Funkcje pomocnicze (SVD, operacje macierzowe)
├── main.cpp               - Program główny z testami i pomiarami
├── Makefile              - Plik kompilacji
└── README.md             - Ten plik
```

## Implementowane algorytmy

### 1. Generator siatki 3D
- Funkcja: `generate3DGridMatrix(k)`
- Generuje macierz o rozmiarze 2^(3k) × 2^(3k)
- Reprezentuje topologię 3D siatki (każdy wierzchołek połączony z 6 sąsiadami)
- Losowe wagi połączeń + dominacja diagonalna dla stabilności

### 2. Kompresja hierarchiczna
- Funkcja: `buildHMatrix(A, maxRank, epsilon)`
- Rekurencyjnie dzieli macierz na 4 ćwiartki
- W liściach: aproksymacja niskiej rangi (SVD)
- Kryterium zatrzymania: mała wartość osobliwa lub mała macierz

### 3. Mnożenie macierz-wektor
- Funkcja: `hMatrixVectorMult(H, x)`
- Algorytm ze slajdu 20
- Liść: Y = U * (V * x)
- Węzeł: rekurencja na synach + scalanie wyników
- Złożoność: O(N log N) zamiast O(N²)

### 4. Dodawanie H-macierzy
- Funkcja: `hMatrixAdd(A, B, maxRank, epsilon)`
- Algorytm ze slajdów 21-22
- Oba liście: konkatenacja U,V + rekompresja
- Oba węzły: rekurencja na odpowiadających synach
- Mixed: podział liścia na strukturę drzewa

### 5. Mnożenie macierz-macierz
- Funkcja: `hMatrixMult(A, B, maxRank, epsilon)`
- Algorytm ze slajdów 23-24
- Oba liście: U_A * (V_A * U_B) * V_B
- Oba węzły: mnożenie blokowe 2×2 z rekursją
- Mixed: podział liścia + rekursja

## Kompilacja i uruchomienie

### Linux/Mac:
```bash
make
./hmatrix
```

### Windows (PowerShell z MinGW):
```powershell
mingw32-make
.\hmatrix.exe
```

### Windows (Visual Studio):
```powershell
g++ -std=c++17 -O3 -o hmatrix.exe main.cpp HMatrix.cpp SupportFunctions.cpp
.\hmatrix.exe
```

## Parametry

Program testuje dla k = 2, 3, 4:
- k=2: macierz 64 × 64
- k=3: macierz 512 × 512
- k=4: macierz 4096 × 4096

Parametry kompresji:
- maxRank = 8
- epsilon = 1e-6

## Wyniki

Program generuje następujące pliki:

1. **hmatrix_structure_k{2,3,4}.txt** - Wizualizacja struktury H-macierzy
2. **timing_vector_mult.txt** - Czasy mnożenia macierz-wektor
3. **timing_matrix_mult.txt** - Czasy mnożenia macierz-macierz
4. **errors.txt** - Błędy aproksymacji

## Analiza złożoności

Program automatycznie oblicza wykładnik β w zależności T = αN^β:
- Dla mnożenia macierz-wektor: oczekiwane β ≈ 1.0-1.5 (liniowa do quasi-liniowej)
- Dla mnożenia macierz-macierz: oczekiwane β ≈ 2.0-2.5 (kwadratowa zamiast sześciennej)

## Struktura węzła H-macierzy

```cpp
struct HNode {
    std::vector<std::shared_ptr<HNode>> sons;  // 4 synów: TL, TR, BL, BR
    int rank;                                   // Ranga aproksymacji
    Matrix U;                                   // Lewy składnik SVD
    Matrix V;                                   // Prawy składnik SVD
    int rows, cols;                            // Wymiary bloku
};
```

## Pseudo-kod kluczowych algorytmów

### Mnożenie macierz-wektor:
```
function hMatrixVectorMult(H, x):
    if H is leaf:
        if H.rank == 0: return zeros
        return H.U * (H.V * x)
    else:
        split x into [x1; x2]
        y1 = mult(H.sons[0], x1) + mult(H.sons[1], x2)
        y2 = mult(H.sons[2], x1) + mult(H.sons[3], x2)
        return [y1; y2]
```

### Mnożenie macierz-macierz:
```
function hMatrixMult(A, B):
    if both leaves:
        return compress(A.U * (A.V * B.U) * B.V)
    else if both internal:
        result.sons[0] = add(mult(A.sons[0], B.sons[0]), mult(A.sons[1], B.sons[2]))
        result.sons[1] = add(mult(A.sons[0], B.sons[1]), mult(A.sons[1], B.sons[3]))
        result.sons[2] = add(mult(A.sons[2], B.sons[0]), mult(A.sons[3], B.sons[2]))
        result.sons[3] = add(mult(A.sons[2], B.sons[1]), mult(A.sons[3], B.sons[3]))
        return result
    else:
        split leaf into 4 blocks and recurse
```

## Autorzy

Marcel Duda, Jan Gawroński

## Data

12.01.2026
