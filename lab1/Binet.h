#pragma once

#include "Mnozenie.h"
#include <memory>

/**
 * Fabryka zwracająca implementację IMnozenie opartą o rekurencyjne mnożenie
 * (Binet / divide and conquer bez pad'owania).
 *
 * Deklaracja tutaj umożliwia użycie implementacji bez ujawniania szczegółów.
 */
std::unique_ptr<IMnozenie> createBinet();