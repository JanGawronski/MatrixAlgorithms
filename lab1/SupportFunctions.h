#pragma once

#include "Mnozenie.h"
#include <iostream>
#include <cstdint>

Matrix createRandomMatrix(int n);
void printSmall(const Matrix& M);

// przeniesione funkcje pomocnicze dla wielu implementacji
Matrix zeroMatrix(int rows, int cols);
Matrix zeroMatrix(int n);
Matrix subMatrix(const Matrix &A, int row, int col, int rows, int cols);
Matrix addMatrix(const Matrix &A, const Matrix &B);
Matrix subtractMatrix(const Matrix &A, const Matrix &B);
Matrix multiplyClassic(const Matrix &A, const Matrix &B);

// Op counter
struct OpCounts {
    std::uint64_t adds = 0;
    std::uint64_t subs = 0;
    std::uint64_t muls = 0;
    std::uint64_t divs = 0;
};

void opCounterReset();
OpCounts opCounterGet();
void opCounterAdd(const OpCounts &c); // optional helper

// Memory accounting (approximate, based on recursive calls)
// current_bytes = currently accounted bytes (sum over active calls of p*r*sizeof(double))
// peak_bytes = maximal observed current_bytes
// active_calls / peak_calls = number of active recursive calls (instantaneous / peak)
struct MemStats {
    std::uint64_t current_bytes = 0;
    std::uint64_t peak_bytes = 0;
    std::uint64_t active_calls = 0;
    std::uint64_t peak_calls = 0;
};

void memCounterReset();
void memCounterEnterCall(std::size_t p, std::size_t r); // account for p*r*sizeof(double)
void memCounterExitCall(std::size_t p, std::size_t r);
MemStats memCounterGet();