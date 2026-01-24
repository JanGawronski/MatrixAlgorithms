@echo off
cd /d c:\Users\Marcel\Documents\studia\MatrixAlgorithms\lab4
g++ -std=c++17 -O3 -o test_vector.exe test_vector.cpp HMatrix.cpp SupportFunctions.cpp
if %ERRORLEVEL% EQU 0 (
    echo.
    echo === Kompilacja zakonczona pomyslnie ===
    echo.
    test_vector.exe
) else (
    echo.
    echo === Blad kompilacji ===
)
pause
