# Test mnożenia macierzy przez wektor [1,2,3,4,5,6,7,8]

## Szybki test (bez H-macierzy)

Otwórz **nowy terminal PowerShell** i uruchom:

```powershell
cd c:\Users\Marcel\Documents\studia\MatrixAlgorithms\lab4
g++ -std=c++17 -o simple_test.exe simple_test.cpp
.\simple_test.exe
```

To pokaże wynik prostego mnożenia macierzy przez wektor bez kompresji.

---

## Test z H-macierzą

Jeśli chcesz zobaczyć działanie H-macierzy:

```powershell
cd c:\Users\Marcel\Documents\studia\MatrixAlgorithms\lab4
g++ -std=c++17 -O3 -o test_vector.exe test_vector.cpp HMatrix.cpp SupportFunctions.cpp
.\test_vector.exe
```

To pokaże:
- Oryginalną macierz
- Wynik mnożenia gęstego (A * x)
- Strukturę H-macierzy (wizualizacja)
- Wynik mnożenia H-macierzy (H * x)
- Błąd aproksymacji

---

## Alternatywnie - użyj pliku batch

Kliknij dwukrotnie plik: `run_test.bat`

---

## Macierz z zadania

Program używa przykładowej macierzy 8×8 o strukturze:

```
[4  1  0  0  1  0  0  0]
[1  4  1  0  0  1  0  0]
[0  1  4  1  0  0  1  0]
[0  0  1  4  0  0  0  1]
[1  0  0  0  4  1  0  0]
[0  1  0  0  1  4  1  0]
[0  0  1  0  0  1  4  1]
[0  0  0  1  0  0  1  4]
```

Wektor: `x = [1, 2, 3, 4, 5, 6, 7, 8]`

---

## Jeśli masz inną macierz z tablicy

Edytuj plik `simple_test.cpp` lub `test_vector.cpp` w linii gdzie zdefiniowana jest macierz A:

```cpp
Matrix A = {
    {4, 1, 0, 0, 1, 0, 0, 0},  // Wiersz 1
    {1, 4, 1, 0, 0, 1, 0, 0},  // Wiersz 2
    // ... itd
};
```

Wpisz wartości z tablicy i ponownie skompiluj.

---

## Oczekiwany wynik dla przykładowej macierzy

Dla macierzy diagonalnie-dominującej z przykładu:

```
A * x ≈ [9, 16, 21, 24, 25, 28, 29, 28]
```

(dokładne wartości zależą od konkretnej macierzy z tablicy)
