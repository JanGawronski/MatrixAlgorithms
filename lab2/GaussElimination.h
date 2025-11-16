#pragma once

#include "Mnozenie.h"

std::pair<Matrix, Matrix> GaussElimination(const Matrix &A, const Matrix &b, std::unique_ptr<IMnozenie> &multImpl);