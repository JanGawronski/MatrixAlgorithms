#pragma once

#include "Mnozenie.h"
#include <memory>

/**
 * Fabryka zwracająca implementację IMnozenie opartą o algorytm AI.
 *
 * Deklaracja tutaj umożliwia użycie implementacji bez ujawniania szczegółów.
 */
std::unique_ptr<IMnozenie> createAI();