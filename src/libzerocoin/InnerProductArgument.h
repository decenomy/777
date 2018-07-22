/**
* @file       InnerProductArgument.h
*
* @brief      InnerProductArgument class for the Zerocoin library.
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

class InnerProductArgument {
public:
    InnerProductArgument(const ZerocoinParams* ZCp): pi(ZKP_MS.size()-1), params(ZCp) {};
    void Prove(const CBigNum y, const CBN_vector a_sets, const CBN_vector b_sets);
    bool Verify(const ZerocoinParams* ZCp, const CBigNum y, CBigNum A, CBigNum B, CBigNum z);
    pair<CBN_vector, CBN_vector> ck_inner_gen(const ZerocoinParams* ZCp, CBigNum y);
    CBN_matrix final_a, final_b;        // final witness
    vector<preChallengeShifts> pi;      // shifted commitments to a_sets and b_sets

private:
    const ZerocoinParams* params;
    CBN_matrix splitIntoSets(const CBN_vector, const int s);
    CBigNum findAkorBk(const CBN_matrix g_sets, const CBN_matrix a_sets, const int k);
    CBigNum findzk(const CBN_matrix a_sets, const CBN_matrix b_sets, const int k);
    void reduction(const int i, CBN_matrix g_sets, CBN_matrix h_sets, CBN_matrix a_sets, CBN_matrix b_sets);
    CBigNum getInnerProductChallenge(const unsigned int i);
    CBN_matrix get_new_gs_hs(const CBN_matrix g_sets, const CBN_vector xPowers, const int m2);
    CBN_matrix get_new_as_bs(const CBN_matrix a_sets, const CBN_vector xPowers, const int m2);
    CBN_matrix final_gs_hs(const CBN_matrix g_sets, const CBN_matrix xPowersList);
    void get_new_A(CBigNum &A2, const CBigNum A, const int i, const CBN_vector xPowersPos, const CBN_vector xPowersNeg);
    void get_new_B(CBigNum &B2, const CBigNum B, const int i, const CBN_vector xPowersPos, const CBN_vector xPowersNeg);
    void get_new_z(CBigNum &z2, const CBigNum z, const int i, const CBN_vector xPowersPos, const CBN_vector xPowersNeg);
    bool testComsCorrect(CBN_matrix ck_sets, const CBigNum C, CBN_matrix sets);
    bool testzCorrect(const CBigNum z);

};

} /* namespace libzerocoin */
