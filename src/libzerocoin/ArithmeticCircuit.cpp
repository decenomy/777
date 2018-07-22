/**
* @file       ArithmeticCircuit.cpp
*
* @brief      ArithmeticCircuit class for the Zerocoin library.
*
* @author     Mary Maller, Jonathan Bootle and Gian Piero Dionisio
* @date       April 2018
*
* @copyright  Copyright 2018 The PIVX Developers
* @license    This project is released under the MIT license.
**/

#include "ArithmeticCircuit.h"

using namespace libzerocoin;

ArithmeticCircuit::ArithmeticCircuit(const ZerocoinParams* p):
            A(ZKP_M, CBN_vector(ZKP_N)),
            B(ZKP_M, CBN_vector(ZKP_N)),
            C(ZKP_M, CBN_vector(ZKP_N)),
            wA(p->ZKP_wA),
            wB(p->ZKP_wB),
            wC(p->ZKP_wC),
            K(p->ZKP_K),
            YPowers(4*ZKP_SERIALSIZE-2),
            YDash(ZKP_N),
            wAj(ZKP_M, CBN_vector(ZKP_N)),
            wBj(ZKP_M, CBN_vector(ZKP_N)),
            wCj(ZKP_M, CBN_vector(ZKP_N)),
            params(p),
            r_bits(ZKP_SERIALSIZE, 0),
            wCoeffA(p->ZKP_wCoeffA),
            wCoeffB(p->ZKP_wCoeffB)
{}

void ArithmeticCircuit::setWireValues(const PrivateCoin& coin)
{
    const CBigNum a = params->coinCommitmentGroup.g;
    const CBigNum b = params->coinCommitmentGroup.h;
    const CBigNum q = params->serialNumberSoKCommitmentGroup.groupOrder;
    /* ---------------------------------- **** WIRE VALUES **** ----------------------------------
     * -------------------------------------------------------------------------------------------
     * Sets wire values (in M*N matrices A, B and C) correctly for a circuit
     * with serial number 'serialNumber' and randomness bits 'r_bits'.
     * Here we have a and b are the group elements used to mint coins.
     */
    serialNumber = coin.getSerialNumber();
    randomness = coin.getRandomness();
    coin.getRandomnessBits(r_bits);

    unsigned int row=0, col=0;
    for(unsigned int i=0; i<ZKP_SERIALSIZE; i++) {
        row = i / ZKP_N;
        col = i % ZKP_N;
        A[row][col] = r_bits[i] % q;
        B[row][col] = (r_bits[i] - CBigNum(1)) % q;
        C[row][col] = CBigNum(0);
    }

    int k=0;
    CBigNum x, product = CBigNum(1);    // efficiency registers
    x = r_bits[0] * (b - CBigNum(1)) + CBigNum(1);
    for(unsigned int i=ZKP_SERIALSIZE; i<2*ZKP_SERIALSIZE-1; i++) {
        row = i / ZKP_N;
        col = i % ZKP_N;
        product = product.mul_mod(x, q);
        A[row][col] = product;
        x = r_bits[k+1] * (b.pow_mod(CBigNum(2).pow(k+1), q) - CBigNum(1)) + CBigNum(1);
        B[row][col] = x % q;
        C[row][col] = A[row][col].mul_mod(B[row][col], q);

        if (i == 2*ZKP_SERIALSIZE-2) {
            A[row][col] = A[row][col].mul_mod(a.pow_mod(serialNumber, q), q);
            C[row][col] = C[row][col].mul_mod(a.pow_mod(serialNumber, q), q);
        }
        k++;
    }
}

void ArithmeticCircuit::setPreConstraints(const ZerocoinParams* params,
        vector<CBN_matrix>& wA, vector<CBN_matrix>& wB, vector<CBN_matrix>& wC, CBN_vector& K,
        CBN_vector& wCoeffA, CBN_vector& wCoeffB)
{
    const CBigNum a = params->coinCommitmentGroup.g;
    const CBigNum b = params->coinCommitmentGroup.h;
    const CBigNum q = params->serialNumberSoKCommitmentGroup.groupOrder;
    /* ---------------------------------- **** CONSTRAINTS **** ----------------------------------
     * -------------------------------------------------------------------------------------------
     * Matrices wA, wB, wC and vector K, specifying constraints that ensure that the circuit
     * is satisfied if and only if Cfinal = a^S b^v.
     *
     */

    wA.resize(4*ZKP_SERIALSIZE-2, CBN_matrix(ZKP_M, CBN_vector(ZKP_N)));
    wB.resize(4*ZKP_SERIALSIZE-2, CBN_matrix(ZKP_M, CBN_vector(ZKP_N)));
    wC.resize(4*ZKP_SERIALSIZE-2, CBN_matrix(ZKP_M, CBN_vector(ZKP_N)));
    K.resize(4*ZKP_SERIALSIZE-2);
    unsigned int i_div_n = 0, i_mod_n = 0;
    unsigned int k_div_n = 0, k_mod_n = 0;
    unsigned int ell_div_n = 0, ell_mod_n = 0;
    unsigned int ell=0;
    int k = 0;
    CBigNum x;             // efficiency reg
    CBN_vector U(ZKP_N);   // Unit vector

    /* Constraints to ensure A[k][l] - B[k][l] = 1 */
    for(unsigned int i=0; i<ZKP_SERIALSIZE; i++) {
        i_div_n = i / ZKP_N;
        i_mod_n = i % ZKP_N;;
        unit_vector(wA[i][i_div_n], i_mod_n);
        unit_vector(U, i_mod_n);
        vectorTimesConstant(wB[i][i_div_n], U, CBigNum(-1), q);
        K[i] = CBigNum(1);
    }

    /* Constraints to ensure C[i] = 0 */
    for(unsigned int i=ZKP_SERIALSIZE; i<2*ZKP_SERIALSIZE; i++) {
        k_div_n = k / ZKP_N;
        k_mod_n = k % ZKP_N;
        unit_vector(wC[i][k_div_n], k_mod_n);
        K[i] = CBigNum(0);
        k++;
    }

    /* Constraints to ensure that B[N+k] = A[k+1] * (b** (2**(k+1)) - 1) +1    */
    k = 1;   ell = ZKP_SERIALSIZE;
    for(unsigned int i=2*ZKP_SERIALSIZE; i<3*ZKP_SERIALSIZE-1; i++) {
        k_div_n = k / ZKP_N;
        k_mod_n = k % ZKP_N;
        ell_div_n = ell / ZKP_N;
        ell_mod_n = ell % ZKP_N;
        unit_vector(U, k_mod_n);
        x = b.pow_mod( CBigNum(2).pow(k), q) - CBigNum(1);
        vectorTimesConstant(wA[i][k_div_n], U, x, q);
        unit_vector(U, ell_mod_n);
        vectorTimesConstant(wB[i][ell_div_n], U, CBigNum(-1), q);
        K[i] = CBigNum(-1) % q;
        k++; ell++;
    }

    /* Constraints to ensure A[N] = A[0] * (b-1) + 1     */
    unit_vector(wA[3*ZKP_SERIALSIZE-1][ZKP_SERIALSIZE/ZKP_N], ZKP_SERIALSIZE%ZKP_N);
    unit_vector(U, 0);
    vectorTimesConstant(wB[3*ZKP_SERIALSIZE-1][0], U, CBigNum(1)-b, q);
    K[3*ZKP_SERIALSIZE-1] = b;

    /* Constraints to ensure A[N + k] = C[N + k -1]     */
    k = ZKP_SERIALSIZE+1; ell = ZKP_SERIALSIZE;
    for(unsigned int i=3*ZKP_SERIALSIZE; i<4*ZKP_SERIALSIZE-3; i++) {
        k_div_n = k/ZKP_N;
        k_mod_n = k%ZKP_N;
        ell_div_n = ell/ZKP_N;
        ell_mod_n = ell%ZKP_N;
        unit_vector(wA[i][k_div_n], k_mod_n);
        unit_vector(U, ell_mod_n);
        vectorTimesConstant(wC[i][ell_div_n], U, CBigNum(-1), q);
        K[i] = 0;
        k++; ell++;
    }

    /* Constraints to ensure A[final] = (a**S) * C[final-1]    */
    k_div_n = (2*ZKP_SERIALSIZE-2) / ZKP_N;
    k_mod_n = (2*ZKP_SERIALSIZE-2) % ZKP_N;
    ell_div_n = (2*ZKP_SERIALSIZE-3) / ZKP_N;
    ell_mod_n = (2*ZKP_SERIALSIZE-3) % ZKP_N;
    unit_vector(wA[4*ZKP_SERIALSIZE-3][k_div_n], k_mod_n);
//    unit_vector(U, ell_mod_n);
//    x = CBigNum(-1).mul_mod(a.pow_mod(serialNumber, q), q);
//    vectorTimesConstant(wC[4*ZKP_SERIALSIZE-3][ell_div_n], U, x, q);
    K[4*ZKP_SERIALSIZE-3] = 0;

    /* Set wCoefficients */
    wCoeffA.resize(4*ZKP_SERIALSIZE-2);
    wCoeffB.resize(4*ZKP_SERIALSIZE-2);
    for(unsigned int i=0; i<wCoeffA.size(); i++) {
        wCoeffA[i] = CBigNum(0);
        wCoeffB[i] = CBigNum(0);
        for(int k=0; k<ZKP_N; k++) {
            wCoeffA[i] = (wCoeffA[i] + wA[i][1][k]) % q;
            wCoeffB[i] = (wCoeffB[i] + wB[i][1][k]) % q;
        }
    }

}

void ArithmeticCircuit::setConstraints(const CBigNum& serialNumber)
{
    const CBigNum q = params->serialNumberSoKCommitmentGroup.groupOrder;
    const CBigNum a = params->coinCommitmentGroup.g;
    CBN_vector U(ZKP_N);   // Unit vector
    unit_vector(U, (2*ZKP_SERIALSIZE-3) % ZKP_N);
    CBigNum x = CBigNum(-1).mul_mod(a.pow_mod(serialNumber, q), q);
    vectorTimesConstant(wC[4*ZKP_SERIALSIZE-3][(2*ZKP_SERIALSIZE-3) / ZKP_N], U, x, q);
}

void ArithmeticCircuit::setYPoly(const CBigNum& y)
{
    /* --------------------------------- **** w-POLYNOMIALS **** ---------------------------------
     * -------------------------------------------------------------------------------------------
     * Matrices wAj, wBj, wCj, specifying  the vector polynomials and array YDash specifying
     * the vector of monomials.
     */
    this->y = y;
    set_YPowers();
    set_YDash();
    set_wABj();
    set_wCj();
    set_Kconst();

}

CBigNum ArithmeticCircuit::sumWiresDotWs(const int i)
{
    CBigNum q = params->serialNumberSoKCommitmentGroup.groupOrder;
    CBigNum sum = CBigNum(0);
    for(unsigned int j=0; j<ZKP_M; j++) {
        sum = (sum + dotProduct(A[j], wA[i][j], q)) % q;
        sum = (sum + dotProduct(B[j], wB[i][j], q)) % q;
        sum = (sum + dotProduct(C[j], wC[i][j], q)) % q;
    }
    return sum;
}

CBigNum ArithmeticCircuit::sumWiresDotWPoly()
{
    CBigNum q = params->serialNumberSoKCommitmentGroup.groupOrder;
    CBigNum x;
    CBigNum sum = CBigNum(0);
    for(unsigned int i=0; i<ZKP_M; i++) {
        x = AiDotBiYDash(i);
        sum = (sum + x.mul_mod(y.pow_mod(i+1,q),q)) % q;
    }
    for(unsigned int i=0; i<ZKP_M; i++) {
        sum = ( sum + dotProduct(A[i], wAj[i], q) ) % q;
        sum = ( sum + dotProduct(B[i], wBj[i], q) ) % q;
        sum = ( sum + dotProduct(C[i], wCj[i], q) ) % q;
    }
    return sum;
}

CBigNum ArithmeticCircuit::AiDotBiYDash(const int i)
{
    CBigNum q = params->serialNumberSoKCommitmentGroup.groupOrder;
    CBigNum res = CBigNum(0);
    for(unsigned int j=0; j<ZKP_N; j++)
        res = (res + A[i][j] * B[i][j] * YDash[j]) % q;
    return res;
}

void ArithmeticCircuit::set_YPowers()
{
    CBigNum q = params->serialNumberSoKCommitmentGroup.groupOrder;
    YPowers[0] = y.pow_mod(CBigNum(ZKP_SERIALSIZE+ZKP_M+1), q);
    for(unsigned int i=1; i<4*ZKP_SERIALSIZE-2; i++)
        YPowers[i] = y.mul_mod(YPowers[i-1], q);
}

void ArithmeticCircuit::set_YDash()
{
    CBigNum q = params->serialNumberSoKCommitmentGroup.groupOrder;
    for(unsigned int i=0; i<ZKP_N; i++)
        YDash[i] = y.pow_mod(ZKP_M*CBigNum(i+1), q);
}

void ArithmeticCircuit::set_wABj()
{
    CBigNum q = params->serialNumberSoKCommitmentGroup.groupOrder;
    CBigNum sumA, sumB;
    for(unsigned int j=0; j<ZKP_M; j++) {
        if (j==1) {
            sumA = CBigNum(0);
            sumB = CBigNum(0);
            for(unsigned int i=0; i<wCoeffA.size(); i++) {
                sumA = (sumA + wCoeffA[i].mul_mod(YPowers[i], q)) % q;
                sumB = (sumB + wCoeffB[i].mul_mod(YPowers[i], q)) % q;
            }
            wAj[j][0] = sumA % q;
            wBj[j][0] = sumB % q;
            fill(wAj[j].begin()+1, wAj[j].end(), CBigNum(0));
            fill(wBj[j].begin()+1, wBj[j].end(), CBigNum(0));
            continue;
        }
        for(unsigned int k=0; k<ZKP_N; k++) {
            sumA = CBigNum(0);
            sumB = CBigNum(0);
            for(unsigned int i=0; i<4*ZKP_SERIALSIZE-2; i++) {
                sumA += wA[i][j][k].mul_mod(YPowers[i], q);
                sumB += wB[i][j][k].mul_mod(YPowers[i], q);
            }
            wAj[j][k] = sumA % q;
            wBj[j][k] = sumB % q;
        }
    }
}


void ArithmeticCircuit::set_wCj()
{
    CBigNum q = params->serialNumberSoKCommitmentGroup.groupOrder;
    CBigNum sum;
    CBigNum y2 = -y.pow_mod(2,q);
    for(unsigned int k=0; k<ZKP_N; k++) {
        sum = - YDash[k].mul_mod(y, q);
        for(unsigned int i=0; i<4*ZKP_SERIALSIZE-2; i++)
            sum += wC[i][0][k].mul_mod(YPowers[i], q);

        wCj[0][k] = sum % q;
        wCj[1][k] = YDash[k].mul_mod(y2, q);
    }
}

void ArithmeticCircuit::set_Kconst()
{
    CBigNum q = params->serialNumberSoKCommitmentGroup.groupOrder;
    Kconst = CBigNum(0);
    for(unsigned int i=0; i<4*ZKP_SERIALSIZE-3; i++)
        Kconst += K[i].mul_mod(YPowers[i], q);
    Kconst = Kconst % q;
}

// verifies correct assignment
void ArithmeticCircuit::check()
{
    const CBigNum a = params->coinCommitmentGroup.g;
    const CBigNum b = params->coinCommitmentGroup.h;
    const CBigNum q = params->serialNumberSoKCommitmentGroup.groupOrder;

    for(unsigned int i=0; i<ZKP_M; i++)
        for(unsigned int j=0; j<ZKP_N; j++)
            if (A[i][j].mul_mod(B[i][j], q) != C[i][j])
                throw std::runtime_error("ArithmeticCircuit::check() error: code 1");
;
    CBigNum logarithm = (a.pow_mod(serialNumber,q)).mul_mod(b.pow_mod(randomness,q),q);
    if (C[ZKP_M-1][0] != logarithm)
        throw std::runtime_error("ArithmeticCircuit::check() error: code 2");

    for(unsigned int i=0; i<4*ZKP_SERIALSIZE-2; i++)
        if(K[i] != sumWiresDotWs(i)) {
            throw std::runtime_error("ArithmeticCircuit::check() error: code 3");
        }

    if (sumWiresDotWPoly() != Kconst)
        throw std::runtime_error("ArithmeticCircuit::check() error: code 4");
}

