# Integracja Lab3 z Lab4

## Podsumowanie

Wszystkie użyteczne komponenty z Lab3 (kompresja macierzy i wizualizacja) zostały zintegrowane z Lab4 (H-macierze).

## Skopiowane pliki z Lab3 do Lab4

### 1. **Compression.h** i **Compression.cpp**
Funkcje kompresji drzewa czwórkowego z SVD:
- `createTree()` - budowa drzewa kompresji
- `reconstructFromTree()` - odtworzenie macierzy z drzewa
- `drawCompression()` - wizualizacja struktury kompresji
- `TreeNode` struktura - kompatybilna z `HNode`
- Funkcje konwersji `treeNodeToHNode()` i `hNodeToTreeNode()`
- **NOWE:** `saveHMatrixVisualizationPNG()` - zapis wizualizacji jako PNG
- **NOWE:** `saveTreeVisualizationPNG()` - zapis wizualizacji TreeNode jako PNG

### 2. **SupportFunctions.cpp** - rozszerzenia
Dodane funkcje z Lab3:
- `createRandomMatrix()` - generowanie losowych macierzy
- `identityMatrix()` - macierz jednostkowa
- Operatory macierzowe: `operator+`, `operator-`, `operator*`
- `negate()` - negacja macierzy
- `compareMatrices()` - porównanie z tolerancją
- `pad()` / `trim()` - dopełnianie/obcinanie
- `printSmall()` - wyświetlanie małych macierzy
- SVD helpers: `vec_dot`, `vec_norm`, `mat_vec_mul`, `transpose_mul`, `power_iteration`
- **NOWE:** `loadImageRGB()` - wczytywanie obrazów RGB
- **NOWE:** `saveImageRGB()` - zapisywanie obrazów RGB
- **NOWE:** `saveImageGray()` - zapisywanie obrazów w skali szarości

### 3. **stb_image.h** i **stb_image_write.h**
Biblioteki STB do obsługi obrazów (PNG, JPG, etc.):
- Umożliwiają wczytywanie i zapisywanie wizualizacji jako plików graficznych
- Używane przez funkcje `loadImageRGB()`, `saveImageRGB()`, `saveImageGray()`
- Single-header libraries - łatwe w integracji

### 4. **HMatrix.h** - zaktualizowany
Dodane deklaracje wszystkich nowych funkcji z Lab3 + funkcje I/O obrazów.

### 5. **Makefile** - zaktualizowany
Uwzględnia kompilację `Compression.cpp` dla wszystkich targetów:
- `hmatrix` - główny program
- `test_vector` - test mnożenia przez wektor
- `simple_test` - prosty test
- **NOWE:** `visualize_example` - przykład wizualizacji

### 6. **visualize_example.cpp** - NOWY
Przykładowy program demonstracyjny pokazujący 3 metody wizualizacji:
- Metoda 1: Bezpośrednia wizualizacja z HNode
- Metoda 2: Konwersja HNode → TreeNode
- Metoda 3: Użycie `createTree()` z Lab3
- Generuje pliki PNG z wizualizacjami struktury H-macierzy

## Struktura kompatybilności

### TreeNode (Lab3) vs HNode (Lab4)

```cpp
// Lab3: TreeNode
struct TreeNode {
    Vector singularValues;
    Matrix U, V;
    TreeNode* topLeft, *topRight, *bottomLeft, *bottomRight;
};

// Lab4: HNode
struct HNode {
    int rank;
    Matrix U, V;
    std::vector<std::shared_ptr<HNode>> sons; // [0]=topLeft, [1]=topRight, [2]=bottomLeft, [3]=bottomRight
};
```

Funkcje konwersji zapewniają pełną kompatybilność między formatami.

## Kompilacja

```bash
# Windows (PowerShell)
g++ -std=c++17 -O3 -o hmatrix.exe main.cpp HMatrix.cpp SupportFunctions.cpp Compression.cpp
g++ -std=c++17 -O3 -o test_vector.exe test_vector.cpp HMatrix.cpp SupportFunctions.cpp Compression.cpp

# Linux/Mac
make all
make test_vector
make simple_test
```

## Użycie wizualizacji z Lab3

```cpp
#include "Compression.h"

// Budowa drzewa kompresji (metoda z Lab3)
TreeNode* tree = createTree(A, rank, epsilon);

// Wizualizacja struktury (do pliku tekstowego)
Matrix vis = drawCompression(tree, rows(A), cols(A));

// Zapis jako PNG (NOWE!)
saveTreeVisualizationPNG(tree, rows(A), cols(A), "structure.png");

// Konwersja do HNode (jeśli potrzebna)
std::shared_ptr<HNode> hnode = treeNodeToHNode(tree);

// Operacje H-macierzowe
Vector result = hMatrixVectorMult(hnode, x);

// Bezpośredni zapis wizualizacji HNode jako PNG
saveHMatrixVisualizationPNG(hnode, "hmatrix_structure.png");
```

### Przykład użycia funkcji obrazowych

```cpp
// Wczytywanie obrazu
auto [R, G, B] = loadImageRGB("input.png");

// Kompresja każdego kanału osobno
TreeNode* tree_R = createTree(R, rank, epsilon);
TreeNode* tree_G = createTree(G, rank, epsilon);
TreeNode* tree_B = createTree(B, rank, epsilon);

// Rekonstrukcja
Matrix R_compressed = reconstructFromTree(tree_R);
Matrix G_compressed = reconstructFromTree(tree_G);
Matrix B_compressed = reconstructFromTree(tree_B);

// Zapis skompresowanego obrazu
saveImageRGB("output_compressed.png", R_compressed, G_compressed, B_compressed);

// Wizualizacja struktury kompresji
Matrix vis_R = drawCompression(tree_R, rows(R), cols(R));
saveImageGray("compression_structure_R.png", vis_R);
```

## Korzyści z integracji

1. **Reużycie kodu** - nie duplikowanie implementacji SVD i kompresji
2. **Spójność** - ta sama metoda kompresji w obu zadaniach
3. **Wizualizacja** - funkcje `drawCompression` działają dla obu struktur
4. **Eksport graficzny** - możliwość zapisu wizualizacji jako PNG
5. **Kompresja obrazów** - możliwość zastosowania H-macierzy do kompresji obrazów
6. **Testowanie** - możliwość porównania wyników Lab3 i Lab4
7. **Rozszerzalność** - łatwe dodawanie nowych operacji macierzowych

## Status

✅ Kompilacja bez błędów
✅ Wszystkie pliki Lab3 zintegrowane
✅ Testy przechodzą (test_vector.exe)
✅ Makefile zaktualizowany
✅ Dokumentacja kompletna
✅ **Wizualizacje PNG działają** (hmatrix_vis_method1.png, method2.png, method3.png)
✅ Biblioteki STB image zintegrowane
✅ Funkcje I/O obrazów dostępne

## Pliki Lab4 po integracji

```
lab4/
├── Compression.cpp          [NOWY - z Lab3 + rozszerzenia PNG]
├── Compression.h            [NOWY - z Lab3 + rozszerzenia PNG]
├── HMatrix.cpp              [zaktualizowany]
├── HMatrix.h                [zaktualizowany - dodane deklaracje + I/O obrazów]
├── SupportFunctions.cpp     [rozszerzony - funkcje z Lab3 + I/O obrazów]
├── stb_image.h              [NOWY - z Lab3]
├── stb_image_write.h        [NOWY - z Lab3]
├── visualize_example.cpp    [NOWY - demo wizualizacji PNG]
├── main.cpp                 [bez zmian]
├── test_vector.cpp          [bez zmian]
├── simple_test.cpp          [bez zmian]
├── Makefile                 [zaktualizowany]
└── INTEGRACJA_LAB3.md       [ten plik]
```

## Wygenerowane pliki wizualizacji

Po uruchomieniu `visualize_example.exe`:
```
├── hmatrix_vis_method1.png  [wizualizacja bezpośrednio z HNode]
├── hmatrix_vis_method2.png  [wizualizacja przez konwersję HNode→TreeNode]
└── hmatrix_vis_method3.png  [wizualizacja przez createTree z Lab3]
```

## Następne kroki

- ✅ Można używać `createTree()` zamiast `buildHMatrix()` dla kompatybilności z Lab3
- ✅ Możliwość porównania wydajności obu metod kompresji
- ✅ Wizualizacja struktury H-macierzy za pomocą funkcji z Lab3
- ✅ Eksport wizualizacji do plików PNG
- ✅ Kompresja obrazów za pomocą H-macierzy
- 📝 Rozszerzenie raport.tex o odniesienia do Lab3 i wizualizacje PNG
- 📝 Dodanie przykładów kompresji obrazów do raportu

## Kompilacja i uruchomienie

```bash
# Windows (PowerShell)
g++ -std=c++17 -O3 -o hmatrix.exe main.cpp HMatrix.cpp SupportFunctions.cpp Compression.cpp
g++ -std=c++17 -O3 -o test_vector.exe test_vector.cpp HMatrix.cpp SupportFunctions.cpp Compression.cpp
g++ -std=c++17 -O3 -o visualize_example.exe visualize_example.cpp HMatrix.cpp SupportFunctions.cpp Compression.cpp

# Uruchomienie przykładu wizualizacji
.\visualize_example.exe

# Linux/Mac
make all
make test_vector
make visualize_example
./visualize_example
```
