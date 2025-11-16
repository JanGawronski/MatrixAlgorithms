#pragma once

#include "Mnozenie.h"

std::pair<Matrix, Matrix> LUfactorization(const Matrix &A, std::unique_ptr<IMnozenie> &multImpl);

double determinantLU(const Matrix &A, std::unique_ptr<IMnozenie> &multImpl);