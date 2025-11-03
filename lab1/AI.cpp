#include "AI.h"
#include "SupportFunctions.h"

#include <stdexcept>
#include <array>
#include <cmath>

Matrix AI::multiply(const Matrix& A, const Matrix& B) {
    int p = static_cast<int>(A.size());                     // rows of A
    int q = p ? static_cast<int>(A[0].size()) : 0;          // cols of A
    int qb = static_cast<int>(B.size());                    // rows of B
    int r = qb ? static_cast<int>(B[0].size()) : 0;         // cols of B
    if (q != qb) throw std::runtime_error("Incompatible dimensions for AI::multiply");


    if (!(p == 4 && q == 5 && r == 5)) {
        throw std::runtime_error("AI::multiply supports only 4x5 * 5x5 matrices");
    }

    // account memory for this 4x5 result
    memCounterEnterCall(4,5,1);

    // prepare h[1..76] (1-based indexing)
    std::array<double, 77> h{};
    OpCounts local_ops{0,0,0,0};

    // --- h1..h75 (same formulas as reference) ---
    // Note: formulas refer only to A rows 0..3 and columns 0..4 (4x5 A) and B 5x5.
    // (Copied from reference; keep ops counting consistent.)
    {////////////////
        double s = - B[1][0] - B[1][4] - B[2][0];
        local_ops.subs += 3;
        h[1] = A[2][1] * s; local_ops.muls += 1;
    }

    {//////////////////
        double sA = A[1][1] + A[1][4] - A[2][4];
        local_ops.adds += 2; local_ops.subs += 1;
        double sB = - B[1][4] - B[4][0];
        local_ops.subs += 2;
        h[2] = sA * sB; local_ops.muls += 1;
    }

    {///////////////////////
        double sA = -A[2][0] - A[3][0] + A[3][1];
        local_ops.subs += 2; local_ops.adds += 1;
        double sB = -B[0][0] + B[1][4];
        local_ops.subs += 1; local_ops.adds += 1;
        h[3] = sA * sB; local_ops.muls += 1;
    }

    {//////////////////
        double sA = A[0][1] + A[0][3] + A[2][3];
        local_ops.adds += 2;
        double sB = -B[1][4] - B[3][0];
        local_ops.subs += 2;
        h[4] = sA * sB; local_ops.muls += 1;
    }

    {//////////////////
        double sA = A[0][4] + A[1][1] + A[1][4];
        local_ops.adds += 2;
        double sB = -B[1][3] + B[4][0];
        local_ops.subs += 1; local_ops.adds += 1;
        h[5] = sA * sB; local_ops.muls += 1;
    }

    {///////////////
        double sA = -A[1][1] - A[1][4] - A[3][4];
        local_ops.subs += 3;
        double sB = B[1][2] + B[4][0];
        local_ops.adds += 1;
        h[6] = sA * sB; local_ops.muls += 1;
    }

    {////////////////
        double sA = -A[0][0] + A[3][0] - A[3][1];
        local_ops.subs += 2; local_ops.adds += 1;
        double sB = B[0][0] + B[1][3];
        local_ops.adds += 1;
        h[7] = sA * sB; local_ops.muls += 1;
    }

    {/////////////////
        double sA = A[2][1] - A[2][2] - A[3][2];
        local_ops.subs += 2;
        double sB = -B[1][2] + B[2][0];
        local_ops.subs += 1; local_ops.adds += 1;
        h[8] = sA * sB; local_ops.muls += 1;
    }

    {///////////
        double sA = -A[0][1] - A[0][3] + A[3][3];
        local_ops.subs += 2; local_ops.adds += 1;
        double sB = B[1][2] + B[3][0];
        local_ops.adds += 1;
        h[9] = sA * sB; local_ops.muls += 1;
    }

    {//////////////////
        double sA = A[1][1] + A[1][4];
        local_ops.adds += 1;
        h[10] = sA * B[4][0]; local_ops.muls += 1;
    }

    {////////////
        double sA = -A[1][0] - A[3][0] + A[3][1];
        local_ops.subs += 2; local_ops.adds += 1;
        double sB = -B[0][0] + B[1][1];
        local_ops.subs += 1; local_ops.adds += 1;
        h[11] = sA * sB; local_ops.muls += 1;
    }

    {////////////
        double sA = A[3][0] - A[3][1];
        local_ops.subs += 1;
        h[12] = sA * B[0][0]; local_ops.muls += 1;
    }

    {////////////////
        double sA = A[0][1] + A[0][3] + A[1][3];
        local_ops.adds += 2;
        double sB = B[1][1] + B[3][0];
        local_ops.adds += 1;
        h[13] = sA * sB; local_ops.muls += 1;
    }

    {//////////////////////
        double sA = A[0][2] - A[2][1] + A[2][2];
        local_ops.subs += 1; local_ops.adds += 1;
        double sB = B[1][3] + B[2][0];
        local_ops.adds += 1;
        h[14] = sA * sB; local_ops.muls += 1;
    }

    {///////////////////
        double sA = -A[0][1] - A[0][3];
        local_ops.subs += 2;
        h[15] = sA * B[3][0]; local_ops.muls += 1;
    }

    {/////////////
        double sA = -A[2][1] + A[2][2];
        local_ops.subs += 1; local_ops.adds += 1;
        h[16] = sA * B[2][0]; local_ops.muls += 1;
    }

    {//////////////
        double sA = A[0][1] + A[0][3] - A[1][0] + A[1][1] - A[1][2] + A[1][3] - A[2][1] + A[2][2] - A[3][0] + A[3][1];
        local_ops.adds += 6; local_ops.subs += 4;
        h[17] = sA * B[1][1]; local_ops.muls += 1;
    }

    {//////////////
        double s = B[0][0] + B[0][1] + B[4][1];
        local_ops.adds += 2;
        h[18] = A[1][0] * s; local_ops.muls += 1;
    }

    {/////////////
        double s = B[2][0] + B[2][1] + B[4][1];
        local_ops.adds += 2;
        h[19] = - A[1][2] * s; local_ops.muls += 1;
    }

    {//////////////////////
        double sA = -A[0][4] + A[1][0] + A[1][2] - A[1][4];
        local_ops.subs += 2; local_ops.adds += 2;
        double sB = -B[0][0] - B[0][1] + B[0][3] - B[4][1];
        local_ops.subs += 3; local_ops.adds += 1;
        h[20] = sA * sB; local_ops.muls += 1;
    }

    {////////////////
        double sA = A[1][0] + A[1][2] - A[1][4];
        local_ops.adds += 1; local_ops.subs += 1;
        h[21] = sA * B[4][1]; local_ops.muls += 1;
    }

    {///////////////////
        double sA = A[0][2] - A[0][3] - A[1][3];
        local_ops.subs += 2;
        double sB = B[0][0] + B[0][1] - B[0][3] - B[2][0] - B[2][1] + B[2][3] + B[3][3];
        local_ops.adds += 3; local_ops.subs += 3;
        h[22] = sA * sB; local_ops.muls += 1;
    }

    {////////////////
        double s = -B[2][0] + B[2][3] + B[3][3];
        local_ops.subs += 1; local_ops.adds += 1;
        h[23] = A[0][2] * s; local_ops.muls += 1;
    }

    {///////////////////////
        double s = -B[3][3] - B[4][0] + B[4][3];
        local_ops.subs += 2; local_ops.adds += 1;
        h[24] = A[0][4] * s; local_ops.muls += 1;
    }
    {/////////////////
        double s = B[0][0] - B[0][3];
        local_ops.subs += 1;
        h[25] = -A[0][0] * s; local_ops.muls += 1;
    }

    {/////////////////////
        double sA = -A[0][2] + A[0][3] + A[0][4];
        local_ops.subs += 1; local_ops.adds += 2;
        h[26] = sA * B[3][3]; local_ops.muls += 1;
    }

    {////////////////
        double sA = A[0][2] - A[2][0] + A[2][2];
        local_ops.subs += 1; local_ops.adds += 1;
        double sB = B[0][0] - B[0][3] + B[0][4] + B[2][4];
        local_ops.subs += 1; local_ops.adds += 2;
        h[27] = sA * sB; local_ops.muls += 1;
    }

    {//////////////////////
        double s = -B[2][4] - B[3][0] - B[3][4];
        local_ops.subs += 2;
        h[28] = - A[2][3] * s; local_ops.muls += 1;
    }

    {////////////////////////
        double s = B[0][0] + B[0][4] + B[2][4];
        local_ops.adds += 2;
        h[29] = A[2][0] * s; local_ops.muls += 1;
    }

    {///////////////
        double sA = A[2][0] - A[2][2] + A[2][3];
        local_ops.subs += 1; local_ops.adds += 1;
        h[30] = sA * B[2][4]; local_ops.muls += 1;
    }

    {///////////////////////
        double sA = -A[0][3] - A[0][4] - A[2][3];
        local_ops.subs += 3;
        double sB = -B[3][3] - B[4][0] + B[4][3] - B[4][4];
        local_ops.subs += 3; local_ops.adds += 1;
        h[31] = sA * sB; local_ops.muls += 1;
    }

    {/////////////////////////
        double sA = A[1][0] + A[3][0] + A[3][3];
        local_ops.adds += 2;
        double sB = B[0][2] - B[3][0] - B[3][1] - B[3][2];
        local_ops.subs += 2;
        h[32] = sA * sB; local_ops.muls += 1;
    }

    {//////////////
        double s = -B[2][0] - B[2][2];
        local_ops.subs += 2;
        h[33] = A[3][2] * s; local_ops.muls += 1;
    }

    {////////////////
        double s = -B[0][2] + B[3][0] + B[3][2];
        local_ops.subs += 1; local_ops.adds += 2;
        h[34] = A[3][3] * s; local_ops.muls += 1;
    }

    {//////////////////
        double s = B[0][2] + B[4][0] + B[4][2];
        local_ops.adds += 2;
        h[35] = -A[3][4] * s; local_ops.muls += 1;
    }

    {/////////////////////
        double sA = A[1][2] - A[1][4] - A[3][4];
        local_ops.subs += 2;
        double sB = B[2][0] + B[2][1] + B[2][2] + B[4][1];
        local_ops.adds += 3;
        h[36] = sA * sB; local_ops.muls += 1;
    }

    {///////////////////
        double sA = -A[3][0] - A[3][3] + A[3][4];
        local_ops.subs += 2; local_ops.adds += 1;
        h[37] = sA * B[0][2]; local_ops.muls += 1;
    }

    {////////////////////
        double sA = -A[1][2] - A[2][0] + A[2][2] - A[2][3];
        local_ops.subs += 3; local_ops.adds += 1;
        double sB = B[2][4] + B[3][0] + B[3][1] + B[3][4];
        local_ops.adds += 3;
        h[38] = sA * sB; local_ops.muls += 1;
    }

    {///////////////////////
        double sA = -A[2][0] - A[3][0] - A[3][3] + A[3][4];
        local_ops.subs += 3; local_ops.adds += 1;
        double sB = B[0][2] + B[4][0] + B[4][2] + B[4][4];
        local_ops.adds += 3;
        h[39] = sA * sB; local_ops.muls += 1;
    }

    {////////////
        double sA = -A[0][2] + A[0][3] + A[0][4] - A[3][3];
        local_ops.subs += 2; local_ops.adds += 2;
        double sB = -B[2][0] - B[2][2] + B[2][3] + B[3][3];
        local_ops.subs += 2; local_ops.adds += 2;
        h[40] = sA * sB; local_ops.muls += 1;
    }

    {//////////////////
        double sA = -A[0][0] + A[3][0] - A[3][4];
        local_ops.subs += 2; local_ops.adds += 1;
        double sB = B[0][2] + B[2][0] + B[2][2] - B[2][3] + B[4][0] + B[4][2] - B[4][3];
        local_ops.adds += 4; local_ops.subs += 2;
        h[41] = sA * sB; local_ops.muls += 1;
    }

    {////////////////////////////
        double sA = -A[1][0] + A[1][4] - A[2][4];
        local_ops.subs += 2; local_ops.adds += 1;
        double sB = -B[0][0] - B[0][1] - B[0][4] + B[3][0] + B[3][1] + B[3][4] - B[4][1];
        local_ops.subs += 4; local_ops.adds += 3;
        h[42] = sA * sB; local_ops.muls += 1;
    }

    {///////////////////
        double s = B[3][0] + B[3][1];
        local_ops.adds += 1;
        h[43] = A[1][3] * s; local_ops.muls += 1;
    }

    {////////////////////////
        double sA = A[1][2] + A[2][1] - A[2][2];
        local_ops.adds += 1; local_ops.subs += 1;
        double sB = B[1][1] - B[2][0];
        local_ops.subs += 1;
        h[44] = sA * sB; local_ops.muls += 1;
    }

    {///////////////////
        double sA = -A[2][2] + A[2][3] - A[3][2];
        local_ops.subs += 2; local_ops.adds += 1;
        double sB = B[2][4] + B[3][0] + B[3][2] + B[3][4] + B[4][0] + B[4][2] + B[4][4];
        local_ops.adds += 6;
        h[45] = sA * sB; local_ops.muls += 1;
    }

    {/////////////////////
        double s = -B[4][0] - B[4][4];
        local_ops.subs += 2;
        h[46] = -A[2][4] * s; local_ops.muls += 1;
    }

    {////////////////////////
        double sA = A[1][0] - A[1][4] - A[2][0] + A[2][4];
        local_ops.adds += 2; local_ops.subs += 2;
        double sB = B[0][0] + B[0][1] + B[0][4] - B[3][0] - B[3][1] - B[3][4];
        local_ops.adds += 3; local_ops.subs += 3;
        h[47] = sA * sB; local_ops.muls += 1;
    }

    {////////////
        double sA = -A[1][2] + A[2][2];
        local_ops.subs += 1; local_ops.adds += 1;
        double sB = B[1][1] + B[2][1] + B[2][4] + B[3][0] + B[3][1] + B[3][4];
        local_ops.adds += 5;
        h[48] = sA * sB; local_ops.muls += 1;
    }

    {/////////////
        double sA = -A[0][0] - A[0][2] + A[0][3] + A[0][4] - A[1][0] - A[1][2] + A[1][3] + A[1][4];
        local_ops.subs += 4; local_ops.adds += 4;
        double sB = -B[0][0] - B[0][1] + B[0][3];
        local_ops.subs += 2; local_ops.adds += 1;
        h[49] = sA * sB; local_ops.muls += 1;
    }

    {///////////
        double sA = -A[0][3] - A[1][3];
        local_ops.subs += 2;
        double sB = B[1][1] - B[2][0] - B[2][1] + B[2][3] - B[3][1] + B[3][3];
        local_ops.adds += 2; local_ops.subs += 3;
        h[50] = sA * sB; local_ops.muls += 1;
    }

    {
        double sB = B[1][0] + B[1][1] - B[4][0];
        local_ops.adds += 1; local_ops.subs += 1;
        h[51] = A[1][1] * sB; local_ops.muls += 1;
        
    }

    {
        double sB = B[0][0] + B[1][0] + B[1][2];
        local_ops.adds += 2;
        h[52] = A[3][1] * sB; local_ops.muls += 1;
    }

    {
        double sB = -B[1][0] + B[1][3] + B[3][0];
        local_ops.subs += 1; local_ops.adds += 2;
        h[53] = -A[0][1] * sB; local_ops.muls += 1;
    }

    {//////////////
        double sA = A[0][1] + A[0][3] - A[1][1] - A[1][4] - A[2][1] + A[2][2] - A[3][1] + A[3][2] - A[3][3] - A[3][4];
        local_ops.adds += 3; local_ops.subs += 7;
        h[54] = sA * B[1][2]; local_ops.muls += 1;
    }

    {
        double sA = A[0][3] - A[3][3];
        local_ops.subs += 1;
        double sB = -B[1][2] + B[2][0] + B[2][2] - B[2][3] + B[3][2] - B[3][3];
        local_ops.adds += 3; local_ops.subs += 3;
        h[55] = sA * sB; local_ops.muls += 1;
    }

    {
        double sA = A[0][0] - A[0][4] - A[3][0] + A[3][4];
        local_ops.adds += 2; local_ops.subs += 2;
        double sB = B[2][0] + B[2][2] - B[2][3] + B[4][0] + B[4][2] - B[4][3];
        local_ops.adds += 3; local_ops.subs += 2;
        h[56] = sA * sB; local_ops.muls += 1;
    }

    {
        double sA = -A[2][0] - A[3][0];
        local_ops.subs += 2;
        double sB = -B[0][2] - B[0][4] - B[1][4] - B[4][0] - B[4][2] - B[4][4];
        local_ops.subs += 6;
        h[57] = sA * sB; local_ops.muls += 1;
    }

    {
        double sA = -A[0][3] - A[0][4] - A[2][3] - A[2][4];
        local_ops.subs += 4;
        double sB = -B[4][0] + B[4][3] - B[4][4];
        local_ops.subs += 2; local_ops.adds += 1;
        h[58] = sA * sB; local_ops.muls += 1;
    }

    {
        double sA = -A[2][2] + A[2][3] - A[3][2] + A[3][3];
        local_ops.subs += 2; local_ops.adds += 2;
        double sB = B[3][0] + B[3][2] + B[3][4] + B[4][0] + B[4][2] + B[4][4];
        local_ops.adds += 5;
        h[59] = sA * sB; local_ops.muls += 1;
    }

    {
        double sA = A[1][4] + A[3][4];
        local_ops.adds += 1;
        double sB = B[1][2] - B[2][0] - B[2][1] - B[2][2] - B[4][1] - B[4][2];
        local_ops.subs += 5;
        h[60] = sA * sB; local_ops.muls += 1;
    }

    {
        double sA = A[0][3] + A[2][3];
        local_ops.adds += 1;
        double sB = B[0][0] - B[0][3] + B[0][4] - B[1][4] - B[3][3] + B[3][4] - B[4][0] + B[4][3] - B[4][4];
        local_ops.adds += 4; local_ops.subs += 5;
        h[61] = sA * sB; local_ops.muls += 1;
    }

    {
        double sA = A[1][0] + A[3][0];
        local_ops.adds += 1;
        double sB = B[0][1] + B[0][2] + B[1][1] - B[3][0] - B[3][1] - B[3][2];
        local_ops.adds += 2; local_ops.subs += 3;
        h[62] = sA * sB; local_ops.muls += 1;
    }

    {
        double sA = -A[2][2] - A[3][2];
        local_ops.subs += 2;
        double sB = -B[1][2] - B[2][2] - B[2][4] - B[3][0] - B[3][2] - B[3][4];
        local_ops.subs += 6;
        h[63] = sA * sB; local_ops.muls += 1;
    }

    {
        double sA = A[0][0] - A[0][2] - A[0][3] + A[2][0] - A[2][2] - A[2][3];
        local_ops.adds += 1; local_ops.subs += 4;
        double sB = B[0][0] - B[0][3] + B[0][4];
        local_ops.subs += 1; local_ops.adds += 1;
        h[64] = sA * sB; local_ops.muls += 1;
    }

    {
        double sA = -A[0][0] + A[3][0];
        local_ops.subs += 1; local_ops.adds += 1;
        double sB = -B[0][2] + B[0][3] + B[1][3] - B[4][0] - B[4][2] + B[4][3];
        local_ops.adds += 3; local_ops.subs += 3;
        h[65] = sA * sB; local_ops.muls += 1;
    }

    {
        double sA = A[0][0] - A[0][1] + A[0][2] - A[0][4] - A[1][1] - A[1][4] - A[2][1] + A[2][2] - A[3][0] + A[3][1];
        local_ops.adds += 3; local_ops.subs += 7;
        h[66] = sA * B[1][3]; local_ops.muls += 1;
    }

    {
        double sA = A[1][4] - A[2][4];
        local_ops.subs += 1;
        double sB = B[0][0] + B[0][1] + B[0][4] - B[1][4] - B[3][0] - B[3][1] - B[3][4] + B[4][1] + B[4][4];
        local_ops.adds += 4; local_ops.subs += 4;
        h[67] = sA * sB; local_ops.muls += 1;
    }

    {
        double sA = A[0][0] + A[0][2] - A[0][3] - A[0][4] -A[3][0] -A[3][2] + A[3][3] + A[3][4];
        local_ops.adds += 4; local_ops.subs += 4;
        double sB = -B[2][0] - B[2][2] + B[2][3];
        local_ops.subs += 2; local_ops.adds += 1;
        h[68] = sA * sB; local_ops.muls += 1;
    }

    {
        double sA = -A[0][2] + A[0][3] - A[1][2] + A[1][3];
        local_ops.subs += 2; local_ops.adds += 2;
        double sB = -B[1][3] - B[2][0] - B[2][1] + B[2][3] - B[4][1] + B[4][3];
        local_ops.subs += 3; local_ops.adds += 2;
        h[69] = sA * sB; local_ops.muls += 1;
    }

    {
        double sA = A[1][2] - A[1][4] + A[3][2] - A[3][4];
        local_ops.adds += 2; local_ops.subs += 2;
        double sB = -B[2][0] - B[2][1] - B[2][2];
        local_ops.subs += 3;
        h[70] = sA * sB; local_ops.muls += 1;
    }

    {
        double sA = -A[2][0] + A[2][2] - A[2][3] + A[2][4] - A[3][0] + A[3][2] - A[3][3] + A[3][4];
        local_ops.adds += 3; local_ops.subs += 5;
        double sB = -B[4][0] - B[4][2] - B[4][4];
        local_ops.subs += 3;
        h[71] = sA * sB; local_ops.muls += 1;
    }

    {
        double sA = -A[1][0] - A[1][3] - A[3][0] - A[3][3];
        local_ops.subs += 4;
        double sB = B[3][0] + B[3][1] + B[3][2];
        local_ops.adds += 2;
        h[72] = sA * sB; local_ops.muls += 1;
    }

    {
        double sA = A[0][2] - A[0][3] - A[0][4] + A[1][2] - A[1][3] - A[1][4];
        local_ops.adds += 4; local_ops.subs += 5;
        double sB = B[0][0] + B[0][1] - B[0][3] + B[1][3] +B[4][1] - B[4][3];
        local_ops.adds += 4; local_ops.subs += 2;
        h[73] = sA * sB; local_ops.muls += 1;
    }

    {
        double sA = A[1][0] - A[1][2] + A[1][3] - A[2][0] + A[2][2] - A[2][3];
        local_ops.adds += 3; local_ops.subs += 3;
        double sB = B[3][0] + B[3][1] + B[3][4];
        local_ops.adds += 2;
        h[74] = sA * sB; local_ops.muls += 1;
    }

    {
        double sA = A[0][1] + A[0][3] - A[1][1] - A[1][4] - A[2][0] + A[2][1] + A[2][3] + A[2][4] - A[3][0] + A[3][1];
        local_ops.adds += 5; local_ops.subs += 5;
        h[75] = - sA * B[1][4]; local_ops.muls += 1;
    }

    {
        double sA = A[0][2] + A[2][2];
        local_ops.adds += 1;
        double sB = -B[0][0] + B[0][3] - B[0][4] + B[1][3] + B[2][3] - B[2][4];
        local_ops.adds += 3; local_ops.subs += 3;
        h[76] = sA * sB; local_ops.muls += 1;
    }

    // --- assemble final 4x5 C from h's ---
    Matrix C = zeroMatrix(4,5);

    {
        double v = -h[10] + h[12] + h[14] - h[15] - h[16] + h[53] + h[5] - h[66] - h[7];
        local_ops.adds += 3; local_ops.subs += 5;
        C[0][0] = v;
    }

    {
        double v = h[10] + h[11] - h[12] + h[13] + h[15] + h[16] - h[17] - h[44] + h[51];
        local_ops.adds += 5; local_ops.subs += 3;
        C[1][0] = v;
    }

    {
        double v = h[10] - h[12] + h[15] + h[16] - h[1] + h[2] + h[3] - h[4] + h[75];
        local_ops.adds += 5; local_ops.subs += 4;
        C[2][0] = v;
    }

    {
        double v = -h[10] + h[12] - h[15] - h[16] + h[52] + h[54] - h[6] - h[8] + h[9];
        local_ops.adds += 3; local_ops.subs += 6;
        C[3][0] = v;
    }

    {
        double v = h[13] + h[15] + h[20] + h[21] - h[22] + h[23] + h[25] - h[43] + h[49] + h[50];
        local_ops.adds += 6; local_ops.subs += 2;
        C[0][1] = v;
    }

    {
        double v = -h[11] + h[12] - h[13] - h[15] - h[16] + h[17] + h[18] - h[19] - h[21] +h[43] + h[44];
        local_ops.adds += 4; local_ops.subs += 7;
        C[1][1] = v;
    }

    {
        double v = -h[16] - h[19] - h[21] - h[28] - h[29] - h[38] + h[42] + h[44] - h[47] + h[48];
        local_ops.adds += 3; local_ops.subs += 7;
        C[2][1] = v;
    }

    {
        double v = h[11] - h[12] - h[18] + h[21] - h[32] + h[33] - h[34] - h[36] + h[62] - h[70];
        local_ops.adds += 4; local_ops.subs += 6;
        C[3][1] = v;
    }

    {
        double v = h[15] + h[23] + h[24] + h[34] - h[37] + h[40] - h[41] + h[55] - h[56] - h[9];
        local_ops.adds += 5; local_ops.subs += 5;
        C[0][2] = v;
    }

    {
        double v = -h[10] + h[19] + h[32] + h[35] + h[36] + h[37] - h[43] - h[60] - h[6] - h[72];
        local_ops.adds += 5; local_ops.subs += 5;
        C[1][2] = v;
    }

    {
        double v = -h[16] - h[28] + h[33] + h[37] - h[39] + h[45] - h[46] + h[63] - h[71] - h[8];
        local_ops.adds += 4; local_ops.subs += 6;
        C[2][2] = v;
    }

    {
        double v = h[10] + h[15] + h[16] - h[33] + h[34] - h[35] - h[37] - h[54] + h[6] + h[8] - h[9];
        local_ops.adds += 5; local_ops.subs += 5;
        C[3][2] = v;
    }

    {
        double v = -h[10] + h[12] + h[14] - h[16] + h[23] + h[24] + h[25] + h[26] +h[5] - h[66] - h[7];
        local_ops.adds += 6; local_ops.subs += 5;
        C[0][3] = v;
    }

    {
        double v = h[10] + h[18] - h[19] + h[20] - h[22] - h[24] - h[26] - h[5] -h[69] + h[73];
        local_ops.adds += 3; local_ops.subs += 7;
        C[1][3] = v;
    }

    {
        double v = -h[14] + h[16] - h[23] - h[26] + h[27] + h[29] + h[31] + h[46] - h[58] + h[76];
        local_ops.adds += 6; local_ops.subs += 4;
        C[2][3] = v;
    }

    {
        double v = h[12] + h[25] + h[26] - h[33] - h[35] - h[40] + h[41] + h[65] - h[68] - h[7];
        local_ops.adds += 5; local_ops.subs += 5;
        C[3][3] = v;
    }

    {
        double v = h[15] + h[24] + h[25] + h[27] - h[28] + h[30] + h[31] - h[4] + h[61] + h[64];
        local_ops.adds += 7; local_ops.subs += 3;
        C[0][4] = v;
    }

    {
        double v = -h[10] - h[18] - h[2] - h[30] - h[38] + h[42] - h[43] + h[46] + h[67] + h[74];
        local_ops.adds += 4; local_ops.subs += 6;
        C[1][4] = v;
    }

    {
        double v = -h[10] + h[12] - h[15] + h[28] + h[29] - h[2] - h[30] -h[3] + h[46] + h[4] - h[75];
        local_ops.adds += 5; local_ops.subs += 7;
        C[2][4] = v;
    }

    {
        double v = -h[12] - h[29] + h[30] - h[34] + h[35] + h[39] + h[3] -h[45] +h[57] + h[59];
        local_ops.adds += 5; local_ops.subs += 5;
        C[3][4] = v;
    }

    // add local operation counts to global counters
    opCounterAdd(local_ops);

    memCounterExitCall(4,5,1);
    return C;
}

std::unique_ptr<IMnozenie> createAI() {
    return std::make_unique<AI>();
}