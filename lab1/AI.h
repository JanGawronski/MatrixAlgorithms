#pragma once

#include "Mnozenie.h"
#include "SupportFunctions.h"
#include <memory>

/**
 * Prosta implementacja IMnozenie nazwana AI.
 * Implementuje klasyczne mnożenie (triple loop) z liczeniem operacji
 * i przybliżonym accountingiem pamięci via memCounterEnter/Exit.
 */
class AI : public IMnozenie {
public:
    Matrix multiply(const Matrix& A, const Matrix& B) override;
    ~AI() override = default;
};

// fabryka
std::unique_ptr<IMnozenie> createAI();