#include "AI.h"
#include "SupportFunctions.h"

#include <algorithm>
#include <stdexcept>
#include <memory>
#include <vector>

namespace {

using Matrix = ::Matrix;

Matrix multiplyRec(const Matrix &A, const Matrix &B) {
    int Arows = static_cast<int>(A.size());
    int Acols = Arows ? static_cast<int>(A[0].size()) : 0;
    int Brows = static_cast<int>(B.size());
    int Bcols = Brows ? static_cast<int>(B[0].size()) : 0;
    if (Acols != 5 || Arows != 4 || Brows != 5 || Bcols != 5) {
        throw std::runtime_error("Works only for 4x5 * 5x5 matrices");
    }

    memCounterEnterCall(static_cast<std::size_t>(Arows), static_cast<std::size_t>(Bcols));

    
    Matrix M = A * B;

    memCounterExitCall(static_cast<std::size_t>(Arows), static_cast<std::size_t>(Bcols));
    return M;
}

} // namespace (internal)


class AIImpl : public IMnozenie {
public:
    Matrix multiply(const Matrix &A, const Matrix &B) override {
        return multiplyRec(A, B);
    }
};

std::unique_ptr<IMnozenie> createAI() {
    return std::make_unique<AIImpl>();
}

