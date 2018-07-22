/**
* @file       ArithmeticCircuit.h
*
* @brief      ArithmeticCircuit class for the Zerocoin library.
*
* @author     Mary Maller, Jonathan Bootle and Gian Piero Dionisio
* @date       April 2018
*
* @copyright  Copyright 2018 The PIVX Developers
* @license    This project is released under the MIT license.
**/

#pragma once
#include "Coin.h"
#include "zkplib.h"

using namespace std;

namespace libzerocoin {

class ArithmeticCircuit {
public:
    ArithmeticCircuit(const ZerocoinParams* p);
    CBN_matrix A;               // left input wires
    CBN_matrix B;               // right input wires
    CBN_matrix C;               // output wires
    vector<CBN_matrix> wA;      // constraints for the left input wires
    vector<CBN_matrix> wB;      // constraints for the right input wires
    vector<CBN_matrix> wC;      // constraints for the output wires
    CBN_vector K;               // constraints vector
    // w-Polynomials
    CBN_vector YPowers;
    CBN_vector YDash;
    CBN_matrix wAj, wBj, wCj;
    CBigNum Kconst;
    void setWireValues(const PrivateCoin& coin);
    static void setPreConstraints(const ZerocoinParams* params, vector<CBN_matrix>& wA, vector<CBN_matrix>& wB, vector<CBN_matrix>& wC, CBN_vector& K, CBN_vector& wCoeffA, CBN_vector& wCoeffB);
    void setConstraints(const CBigNum& serialNumber);
    void setYPoly(const CBigNum& y);
    const CBigNum& getSerialNumber() const { return this->serialNumber; }
    const CBigNum& getRandomness() const { return this->randomness; }
    CBigNum sumWiresDotWs(const int i);    // Evaluate the sums in Equation (1) of the paper
    CBigNum sumWiresDotWPoly();            // Evaluate the sums in Equation (2) of the paper
    CBigNum AiDotBiYDash(const int i);     // Evaluate dotProduct(A[i], hadamard(B[i], YDash)
    void check();                          // perform tests on circuit assignment
private:
    const ZerocoinParams* params;
    CBigNum serialNumber;       // coin serial number S
    CBigNum randomness;         // coin randomness v
    vector<int> r_bits;               // randomness binary decomposition
    CBigNum y;                        // indeterminate of the polynomial equation
    CBN_vector wCoeffA, wCoeffB;
    void set_YPowers();
    void set_YDash();
    void w1PolynomialCoefficients();
    void set_wABj();
    void set_wCj();
    void set_Kconst();
};

} /* namespace libzerocoin */
