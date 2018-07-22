/**
* @file       InnerProductArgument.cpp
*
* @brief      InnerProductArgument class for the Zerocoin library.
*
* @author     Mary Maller, Jonathan Bootle and Gian Piero Dionisio
* @date       April 2018
*
* @copyright  Copyright 2018 The PIVX Developers
* @license    This project is released under the MIT license.
**/
#include "hash.h"
#include "InnerProductArgument.h"
//#include <time.h>

using namespace libzerocoin;

void InnerProductArgument::Prove(const CBigNum y, const CBN_vector a_sets, const CBN_vector b_sets)
{
    // ----------------------------- **** INNER-PRODUCT PROVE **** ------------------------------
    // ------------------------------------------------------------------------------------------
    // Algorithm used by the signer to prove the inner product argument for two vector polynomials
    // @param   y                       :  value used to create commitment keys ck_inner_{g,h}
    // @param   a_sets, b_sets          :  The two length N+PADS vectors that are committed to
    // @init    final_a, final_b        :  final witness
    // @init    pi                      :  final shifted commitments
/*
    clock_t step_start_time, step_total_time;
    clock_t start_time, total_time;
    cout << "INNER PRODUCT PROVE" << endl;
*/
    const CBigNum q = params->serialNumberSoKCommitmentGroup.groupOrder;
    const CBigNum p = params->serialNumberSoKCommitmentGroup.modulus;

    pair<CBN_vector, CBN_vector> resultSets = ck_inner_gen(params, y);
    CBN_vector ck_inner_g = resultSets.first;
    CBN_vector ck_inner_h = resultSets.second;

    int mu = ZKP_MS.size();
    int s = ZKP_MS[mu-1];

    // [s][M/s] matrices
    CBN_matrix g_sets = splitIntoSets(ck_inner_g, s);
    CBN_matrix h_sets = splitIntoSets(ck_inner_h, s);
    final_a = splitIntoSets(a_sets, s);
    final_b = splitIntoSets(b_sets, s);

    CBN_vector xPowersPos(s+1), xPowersNeg(s+1);
    xPowersPos[0] = xPowersNeg[0] = CBigNum(1);

    CBigNum x;
    unsigned int i = 0;
    while(mu > 1){
 /*
        cout << "- iteration n. " << i << endl;
        cout << " (sets dim: " << g_sets.size() << "*" << g_sets[0].size() << ")" << endl;
        step_start_time = clock();
        start_time = clock();
*/
        reduction(i, g_sets, h_sets, final_a, final_b);
/*
        total_time = clock() - start_time;
        cout << "* reduction completed in " << total_time*1.0/CLOCKS_PER_SEC << " sec" << endl;
        start_time = clock();
*/
        x = getInnerProductChallenge(i);
/*
        total_time = clock() - start_time;
        cout << "* got x challenge in " << total_time*1.0/CLOCKS_PER_SEC << " sec" << endl;
        start_time = clock();
*/
        xPowersPos[1] = x;
        xPowersNeg[1] = x.pow_mod(-1,q);
        for(int j=2; j<=ZKP_MS[mu-2]; j++) {
            xPowersPos[j] = xPowersPos[j-1].mul_mod(x, q);
            xPowersNeg[j] = xPowersNeg[j-1].mul_mod(xPowersNeg[1], q);
        }
/*
        total_time = clock() - start_time;
        cout << "* precomputed x powers in " << total_time*1.0/CLOCKS_PER_SEC << " sec" << endl;
        start_time = clock();
*/
        g_sets = get_new_gs_hs(g_sets, xPowersNeg, ZKP_MS[mu-2]);
        h_sets = get_new_gs_hs(h_sets, xPowersPos, ZKP_MS[mu-2]);
/*
        total_time = clock() - start_time;
        cout << "* got new g_sets and h_sets in " << total_time*1.0/CLOCKS_PER_SEC << " sec" << endl;
        start_time = clock();
*/
        final_a = get_new_as_bs(final_a, xPowersPos, ZKP_MS[mu-2]);
        final_b = get_new_as_bs(final_b, xPowersNeg, ZKP_MS[mu-2]);
/*
        total_time = clock() - start_time;
        cout << "* got new a_sets and b_sets in " << total_time*1.0/CLOCKS_PER_SEC << " sec" << endl;

        step_total_time = clock() - step_start_time;
        cout << "---> total iteration time: " << step_total_time*1.0/CLOCKS_PER_SEC << " sec" << endl << endl;
*/
        mu -= 1;
        i++;
    }

}


bool InnerProductArgument::Verify(const ZerocoinParams* ZCp, const CBigNum y, CBigNum A, CBigNum B, CBigNum z)
{
    // ---------------------------- **** INNER-PRODUCT VERIFY **** ------------------------------
    // ------------------------------------------------------------------------------------------
    // Algorithm used to verify the inner product proof
    // @param   y       :  value used to create commitment keys ck_inner_{g,h}
    // @param   A, B    :  commitments (to a_sets and b_sets)
    // @param   z       :  target value
    // @return  bool    :  result of the verification
/*
    clock_t step_start_time, step_total_time;
    clock_t start_time, total_time;
    cout << "INNER PRODUCT VERIFY" << endl;
*/
    const CBigNum q = params->serialNumberSoKCommitmentGroup.groupOrder;
    const CBigNum p = params->serialNumberSoKCommitmentGroup.modulus;

    pair<CBN_vector, CBN_vector> resultSets = ck_inner_gen(params, y);
    CBN_vector ck_inner_g = resultSets.first;
    CBN_vector ck_inner_h = resultSets.second;

    int mu = ZKP_MS.size();
    int s = ZKP_MS[mu-1];

    // [s][M/s] matrices
    CBN_matrix g_sets = splitIntoSets(ck_inner_g, s);
    CBN_matrix h_sets = splitIntoSets(ck_inner_h, s);

    CBigNum x;
    CBigNum A2, B2, z2;

    CBN_vector xPowersPos(s+1), xPowersNeg(s+1);
    xPowersPos[0] = xPowersNeg[0] = CBigNum(1);

    CBN_matrix xPowersPositiveList, xPowersNegativeList;

    unsigned int i = 0;
    while(mu > 1) {
/*
        cout << "- iteration n. " << i << endl;
        cout << " (sets dim: " << g_sets.size() << "*" << g_sets[0].size() << ")" << endl;
        step_start_time = clock();
        start_time = clock();
*/
        x = getInnerProductChallenge(i);
/*
        total_time = clock() - start_time;
        cout << "* got x challenge in " << total_time*1.0/CLOCKS_PER_SEC << " sec" << endl;
        start_time = clock();
*/
        xPowersPos[1] = x;
        xPowersNeg[1] = x.pow_mod(-1,q);
        for(int j=2; j<=ZKP_MS[mu-2]; j++) {
            xPowersPos[j] = xPowersPos[j-1].mul_mod(x, q);
            xPowersNeg[j] = xPowersNeg[j-1].mul_mod(xPowersNeg[1], q);
        }
/*
        total_time = clock() - start_time;
        cout << "* precomputed x powers in " << total_time*1.0/CLOCKS_PER_SEC << " sec" << endl;
        start_time = clock();
*/
        //g_sets = get_new_gs_hs(g_sets, xPowersNeg, ZKP_MS[mu-2]);
        //h_sets = get_new_gs_hs(h_sets, xPowersPos, ZKP_MS[mu-2]);
        xPowersPositiveList.push_back(CBN_vector(xPowersPos.begin()+1, xPowersPos.end()));
        xPowersNegativeList.push_back(CBN_vector(xPowersNeg.begin()+1, xPowersNeg.end()));
/*
        total_time = clock() - start_time;
        cout << "* got new g_sets and h_sets in " << total_time*1.0/CLOCKS_PER_SEC << " sec" << endl;
        start_time = clock();
*/
        get_new_A(A2, A, i, xPowersPos, xPowersNeg);
        get_new_B(B2, B, i, xPowersPos, xPowersNeg);
        get_new_z(z2, z, i, xPowersPos, xPowersNeg);
/*
        total_time = clock() - start_time;
        cout << "* got new A, B and z in " << total_time*1.0/CLOCKS_PER_SEC << " sec" << endl;

        step_total_time = clock() - step_start_time;
        cout << "---> total iteration time: " << step_total_time*1.0/CLOCKS_PER_SEC << " sec" << endl << endl;
*/
        A = A2;
        B = B2;
        z = z2;
        mu--;
        i++;
    }

    g_sets = final_gs_hs(g_sets, xPowersNegativeList);
    h_sets = final_gs_hs(h_sets, xPowersPositiveList);

    bool b = testComsCorrect(g_sets, A, final_a);
    b = b && testComsCorrect(h_sets, B, final_b);
    b = b && testzCorrect(z);
    return b;
}


// Initialize sets for inner product
pair<CBN_vector, CBN_vector> InnerProductArgument::ck_inner_gen(const ZerocoinParams* ZCp, CBigNum y)
{
    const IntegerGroupParams* SoKgroup = &(ZCp->serialNumberSoKCommitmentGroup);
    const CBigNum q = SoKgroup->groupOrder;
    const CBigNum p = SoKgroup->modulus;
    CBN_vector ck_inner_g;
    CBN_vector ck_inner_h;

    CBigNum exp = CBigNum(1);
    CBigNum ym = y.pow_mod(-ZKP_M, q);

    for(int j=0; j<(ZKP_N+ZKP_PADS); j++) {
        ck_inner_g.push_back(SoKgroup->gis[j]);
        exp = exp.mul_mod(ym,q);
        ck_inner_h.push_back(SoKgroup->gis[j].pow_mod(exp,p));
    }

    return make_pair(ck_inner_g, ck_inner_h);
}


// Splits an [N1]-vector v into an [s][N1/s]-matrix g_sets
CBN_matrix InnerProductArgument::splitIntoSets(const CBN_vector v, const int s)
{
    const int N1 = v.size();
    const int N2 = N1/s;

    // assert s divides M
    if(N1 % s != 0) throw std::runtime_error("error: s does not divide M");

    // allocate g_sets
    CBN_matrix g_sets(s, vector<CBigNum>(N2));

    for(int i=0; i<s; i++) {
        for(int j=0; j<N2; j++)
            g_sets[i][j] = v[i*N2 + j];
    }

    return g_sets;
}


// Functions for reduction when mu > 1
CBigNum InnerProductArgument::findAkorBk(const CBN_matrix g_sets, const CBN_matrix a_sets, const int k)
{
    const CBigNum p = params->serialNumberSoKCommitmentGroup.modulus;
    const int M1 = g_sets.size();
    const int N1 = g_sets[0].size();
    const int r1 = max(0,-k);
    const int r2 = min(M1,M1-k);
    CBigNum Ak = CBigNum(1);

    for(int i=r1; i<r2; i++)
        for(int j=0; j<N1; j++)
        Ak = Ak.mul_mod(g_sets[i][j].pow_mod(a_sets[i+k][j],p),p);

    return Ak;
}

CBigNum InnerProductArgument::findzk(const CBN_matrix a_sets, const CBN_matrix b_sets, const int k)
{
    const CBigNum q = params->serialNumberSoKCommitmentGroup.groupOrder;
    const int M1 = a_sets.size();
    const int r1 = max(0,-k);
    const int r2 = min(M1,M1-k);
    CBigNum zk = CBigNum(0);

    for(int i=r1; i<r2; i++)
        zk = ( zk + dotProduct(a_sets[i], b_sets[i+k], q) ) % q;

    return zk;
}

// Reduction used in innerProductProve when mu > 1
void InnerProductArgument::reduction(const int i, CBN_matrix g_sets, CBN_matrix h_sets, CBN_matrix a_sets, CBN_matrix b_sets)
{
    preChallengeShifts pcs;
    const int M = a_sets.size();

    for(int k=1; k<M; k++) {
        pcs.nAks.push_back(findAkorBk(g_sets, a_sets, -k));
        pcs.pAks.push_back(findAkorBk(g_sets, a_sets, k));

        pcs.nBks.push_back(findAkorBk(h_sets, b_sets, -k));
        pcs.pBks.push_back(findAkorBk(h_sets, b_sets, k));

        pcs.nzks.push_back(findzk(a_sets, b_sets, -k));
        pcs.pzks.push_back(findzk(a_sets, b_sets, k));
    }

    pi[i] = pcs;
}

// Challenge required by both prover and verifier in inner product argument.
CBigNum InnerProductArgument::getInnerProductChallenge(const unsigned int i)
{
    const CBigNum q = params->serialNumberSoKCommitmentGroup.groupOrder;
    const unsigned int M = pi[i].nAks.size();
    CHashWriter1024 hasher(0,0);
    for(unsigned int j=0; j<M; j++) {
        hasher << pi[i].nAks[j].ToString() << pi[i].pAks[j].ToString();
        hasher << pi[i].nBks[j].ToString() << pi[i].pBks[j].ToString();
        hasher << pi[i].nzks[j].ToString() << pi[i].pzks[j].ToString();
    }

    return CBigNum(hasher.GetHash()) % q;
}

// Find the new gs and hs in the inner product reduction
CBN_matrix InnerProductArgument::get_new_gs_hs(const CBN_matrix g_sets,
        const CBN_vector xPowers, const int m2)
{
    const CBigNum q = params->serialNumberSoKCommitmentGroup.groupOrder;
    const CBigNum p = params->serialNumberSoKCommitmentGroup.modulus;
    const int M1 = g_sets.size();
    const int N1 = g_sets[0].size();
    CBN_vector new_gs(N1);
    CBigNum new_g;

    for(int j=0; j<N1; j++) {
        new_g = CBigNum(1);
        for(int i=0; i<M1; i++)
            new_g = new_g.mul_mod(g_sets[i][j].pow_mod(xPowers[i+1],p),p);
        new_gs[j] = new_g;
    }

    return splitIntoSets(new_gs, m2);
}


// Find the new gs and hs in the inner product reduction (verify)
CBN_matrix InnerProductArgument::final_gs_hs(const CBN_matrix g_sets, const CBN_matrix xPowersList)
{
    const CBigNum q = params->serialNumberSoKCommitmentGroup.groupOrder;
    const CBigNum p = params->serialNumberSoKCommitmentGroup.modulus;

    CBN_matrix X(xPowersList.rbegin(), xPowersList.rend());

    const int N1 = g_sets[0].size();
    CBigNum g0 = CBigNum(1);
    CBigNum g1 = CBigNum(1);

    CBN_matrix exponents0(2, CBN_vector(2, CBigNum(0)));
    CBN_matrix exponents1(2, CBN_vector(2, CBigNum(0)));
    CBN_matrix exponents2(2, CBN_vector(2, CBigNum(0)));
    CBN_matrix exponents3(1, CBN_vector(2, CBigNum(0)));

    exponents0[0][0] = X[1][0].mul_mod(X[0][0], q);
    exponents0[0][1] = X[1][0].mul_mod(X[0][1], q);
    exponents0[1][0] = X[1][1].mul_mod(X[0][0], q);
    exponents0[1][1] = X[1][1].mul_mod(X[0][1], q);
    exponents1[0][0] = X[3][0].mul_mod(X[2][0], q);
    exponents1[0][1] = X[3][0].mul_mod(X[2][1], q);
    exponents1[1][0] = X[3][1].mul_mod(X[2][0], q);
    exponents1[1][1] = X[3][1].mul_mod(X[2][1], q);
    exponents2[0][0] = X[5][0].mul_mod(X[4][0], q);
    exponents2[0][1] = X[5][0].mul_mod(X[4][1], q);
    exponents2[1][0] = X[5][1].mul_mod(X[4][0], q);
    exponents2[1][1] = X[5][1].mul_mod(X[4][1], q);
    exponents3[0][0] = X[6][0];
    exponents3[0][1] = X[6][1];


    vector<int> binary_i;
    CBigNum expo, expo0, expo1;

    for(int i=0; i<N1/2; i++) {
        binary_lookup(binary_i, i);

        expo = exponents0[binary_i[1]][binary_i[0]].mul_mod(
                exponents1[binary_i[3]][binary_i[2]].mul_mod(
                        exponents2[binary_i[5]][binary_i[4]].mul_mod(
                                X[6][binary_i[6]],q),q),q);

        expo0 = (X[ZKP_MS.size()-2][0].mul_mod(expo, q));
        expo1 = (X[ZKP_MS.size()-2][1].mul_mod(expo, q));

        g0 = g0.mul_mod(g_sets[0][2*i].pow_mod(expo0, p), p);
        g0 = g0.mul_mod(g_sets[1][2*i].pow_mod(expo1, p), p);
        g1 = g1.mul_mod(g_sets[0][2*i+1].pow_mod(expo0, p), p);
        g1 = g1.mul_mod(g_sets[1][2*i+1].pow_mod(expo1, p), p);
    }

    CBN_matrix final_g_sets(2, CBN_vector(1, CBigNum(0)));

    final_g_sets[0][0] = g0;
    final_g_sets[1][0] = g1;

    return final_g_sets;
}


// Get new as and bs in inner product reduction
CBN_matrix InnerProductArgument::get_new_as_bs(const CBN_matrix a_sets,
        const CBN_vector xPowers, const int m2)
{
    const CBigNum q = params->serialNumberSoKCommitmentGroup.groupOrder;
    const int M1 = a_sets.size();
    const int N1 = a_sets[0].size();
    CBN_vector a2(N1);
    CBigNum aj;

    for(int j=0; j<N1; j++) {
        aj = CBigNum(0);
        for(int i=0; i<M1; i++) {
            aj += a_sets[i][j].mul_mod(xPowers[i+1],q);
            aj %= q;
        }
        a2[j] = aj;
    }

    return splitIntoSets(a2, m2);
}


// Get new commitments A,B in inner product argument
void InnerProductArgument::get_new_A(CBigNum &A2, const CBigNum A, const int i,
        const CBN_vector xPowersPos, const CBN_vector xPowersNeg)
{
    const CBigNum q = params->serialNumberSoKCommitmentGroup.groupOrder;
    const CBigNum p = params->serialNumberSoKCommitmentGroup.modulus;
    const int M1 = pi[i].nAks.size();
    A2 = A;

    for(int j=0; j<M1; j++) {
        A2 = A2.mul_mod(pi[i].nAks[j].pow_mod(xPowersNeg[j+1],p),p);
        A2 = A2.mul_mod(pi[i].pAks[j].pow_mod(xPowersPos[j+1],p),p);
    }
}

void InnerProductArgument::get_new_B(CBigNum &B2, const CBigNum B, const int i,
        const CBN_vector xPowersPos, const CBN_vector xPowersNeg)
{
    const CBigNum q = params->serialNumberSoKCommitmentGroup.groupOrder;
    const CBigNum p = params->serialNumberSoKCommitmentGroup.modulus;
    const int M1 = pi[i].nBks.size();
    B2 = B;

    for(int j=0; j<M1; j++) {
        B2 = B2.mul_mod(pi[i].nBks[j].pow_mod(xPowersPos[j+1],p),p);
        B2 = B2.mul_mod(pi[i].pBks[j].pow_mod(xPowersNeg[j+1],p),p);
    }
}

void InnerProductArgument::get_new_z(CBigNum &z2, const CBigNum z, const int i,
        const CBN_vector xPowersPos, const CBN_vector xPowersNeg)
{
    const CBigNum q = params->serialNumberSoKCommitmentGroup.groupOrder;
    const int M1 = pi[i].nzks.size();
    z2 = z;

    for(int j=0; j<M1; j++) {
        z2 += pi[i].nzks[j].mul_mod(xPowersPos[j+1],q);
        z2 += pi[i].pzks[j].mul_mod(xPowersNeg[j+1],q);
        z2 %= q;
    }

}

// check commitment C under commitment key ck_set opens to sets (final_a or final_b)
bool InnerProductArgument::testComsCorrect(CBN_matrix ck_sets, const CBigNum C, CBN_matrix sets)
{
    const CBigNum p = params->serialNumberSoKCommitmentGroup.modulus;
    int M1 = sets.size();
    int N1 = sets[0].size();

    CBigNum Ctest = CBigNum(1);
    for(int i=0; i<M1; i++) for(int j=0; j<N1; j++)
        Ctest = Ctest.mul_mod(ck_sets[i][j].pow_mod(sets[i][j],p),p);

    if(Ctest != C) {
        LogPrintf("InnerProductArgument::Verify() - commitment check failed\nC = %s\nCtest = %s", C.ToString(), Ctest.ToString());

        return false;
    }

    return true;
}

// check z is inner product final_a and final_b
bool InnerProductArgument::testzCorrect(const CBigNum z)
{
    const CBigNum q = params->serialNumberSoKCommitmentGroup.groupOrder;
    int M1 = final_a.size();

    CBigNum ztest = CBigNum(0);
    for(int i=0; i<M1; i++) {
        ztest += dotProduct(final_a[i], final_b[i], q);
        ztest %= q;
    }

    if(ztest != z) {
        LogPrintf("InnerProductArgument::Verify() - z check failed\nz = %s\nztest = %s", z.ToString(), ztest.ToString());
        return false;
    }

    return true;
}
