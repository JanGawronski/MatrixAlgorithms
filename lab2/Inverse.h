#pragma once

#include "Mnozenie.h"

Matrix inverse(const Matrix &A, std::unique_ptr<IMnozenie> &multImpl);