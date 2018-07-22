/**
* @file       SerialNumberSoK_small.cpp
*
* @brief      SerialNumberSoK_small class for the Zerocoin library.
*
* @author     Mary Maller, Jonathan Bootle and Gian Piero Dionisio
* @date       April 2018
*
* @copyright  Copyright 2018 The PIVX Developers
* @license    This project is released under the MIT license.
**/
#include <streams.h>
#include "ArithmeticCircuit.h"
#include "SerialNumberSoK_small.h"
//#include <time.h>

namespace libzerocoin {

SerialNumberSoK_small::SerialNumberSoK_small(const ZerocoinParams* ZCp) :
                params(ZCp),
                ComA(ZKP_M),
                ComB(ZKP_M),
                ComC(ZKP_M),
                polyComm(ZCp),
                innerProduct(ZCp)
{ }


SerialNumberSoK_small::SerialNumberSoK_small(const ZerocoinParams* ZCp, const PrivateCoin& coin,
        const Commitment& commitmentToCoin, uint256 msghash) :
                params(ZCp),
                ComA(ZKP_M),
                ComB(ZKP_M),
                ComC(ZKP_M),
                polyComm(ZCp),
                innerProduct(ZCp)

{
    // ---------------------------------- **** SoK PROVE **** -----------------------------------
    // ------------------------------------------------------------------------------------------
    // Specifies how a spender should produce the signature of knowledge on a message msghash that
    // he knows v such that commitmentToCoin is a commitment to a^S b^v.
    // @param   coin                :  The PrivateCoin ww are committing to
    // @param   commitmentToCoin    :  commitment (y1)
    // @param   msghash             :  a message hash to sign
    // @init    SoK

    const CBigNum a = params->coinCommitmentGroup.g;
    const CBigNum b = params->coinCommitmentGroup.h;
    const CBigNum g = params->serialNumberSoKCommitmentGroup.g;
    const CBigNum h = params->serialNumberSoKCommitmentGroup.h;
    const CBigNum q = params->serialNumberSoKCommitmentGroup.groupOrder;
    const CBigNum p = params->serialNumberSoKCommitmentGroup.modulus;
    const CBigNum y1 = commitmentToCoin.getCommitmentValue();
    const CBigNum S = coin.getSerialNumber();
    const CBigNum v = coin.getRandomness();
    const int m =  ZKP_M;
    const int n =  ZKP_N;
    const int m1dash =  ZKP_M1DASH;
    const int m2dash =  ZKP_M2DASH;
    const int ndash =  ZKP_NDASH;
    const int pads =  ZKP_PADS;

    // ****************************************************************************
    // ********************** STEP 1: Generate Commitments ************************
    // ****************************************************************************

    // Select blinding vectors alpha, beta, gamma, D, delta
    CBN_vector f_alpha(m);
    CBN_vector f_beta(m);
    CBN_vector f_gamma(m);
    CBN_vector D(n);
    CBigNum f_delta = CBigNum::randBignum(q);

    random_vector_mod(f_alpha, q);
    random_vector_mod(f_beta, q);
    random_vector_mod(f_gamma, q);
    random_vector_mod(D, q);

    // set arithmetic circuit wire values and constraints
    ArithmeticCircuit circuit(params);
    circuit.setWireValues(coin);
    circuit.setConstraints(coin.getSerialNumber());

    // Commit to the assignment of the circuit
    for(unsigned int i=0; i<m; i++) {
        ComA[i] = pedersenCommitment(params, circuit.A[i], f_alpha[i]);
        ComB[i] = pedersenCommitment(params, circuit.B[i], f_beta[i]);
        ComC[i] = pedersenCommitment(params, circuit.C[i], f_gamma[i]);
    }
    ComD = pedersenCommitment(params, D, f_delta);

    // replace commitment y1 and blind value r
    ComC[m-1] = y1;
    f_gamma[m-1] = commitmentToCoin.getRandomness();


    // ****************************************************************************
    // ************* STEP 2: Challenge component + eval w-polynomials *************
    // ****************************************************************************

    CHashWriter1024 hasher(0,0);
    hasher << msghash << ComD.ToString();

    for(unsigned int i=0; i<m; i++)
        hasher << ComA[i].ToString() << ComB[i].ToString() << ComC[i].ToString();

    // get the challenge component y
    CBigNum y = CBigNum(hasher.GetHash() )% q;

    // set circuit w-Polynomials
    circuit.setYPoly(y);

    // verify correct assignment of circuit values
    // !TODO: skip this for efficiency?
    //circuit.check();


    // ****************************************************************************
    // ************************ STEP 3: Laurent polynomial ************************
    // ****************************************************************************

    // rPoly and sPoly
    CBN_matrix rPolyPositive(2*m+2, CBN_vector(n)), rPolyNegative(2*m+2, CBN_vector(n));
    CBN_matrix sPolyPositive(2*m+2, CBN_vector(n)), sPolyNegative(2*m+2, CBN_vector(n));

    fill(rPolyPositive[0].begin(), rPolyPositive[0].end(), CBigNum(0));
    fill(rPolyNegative[0].begin(), rPolyNegative[0].end(), CBigNum(0));
    fill(sPolyPositive[0].begin(), sPolyPositive[0].end(), CBigNum(0));
    fill(sPolyNegative[0].begin(), sPolyNegative[0].end(), CBigNum(0));

    for(int i=1; i<m+1; i++) {
        rPolyNegative[i] = circuit.B[i-1];
        sPolyPositive[i] = circuit.wBj[i-1];

        for(int j=0; j<n; j++) {
            rPolyPositive[i][j] = circuit.A[i-1][j].mul_mod(y.pow_mod(i,q),q);
            sPolyNegative[i][j] = circuit.wAj[i-1][j].mul_mod(y.pow_mod(-i,q),q);
        }
    }

    for(int i=m+1; i<2*m+1; i++) {
        rPolyPositive[i] = circuit.C[i-1-m];
        fill(rPolyNegative[i].begin(), rPolyNegative[i].end(), CBigNum(0));
        fill(sPolyPositive[i].begin(), sPolyPositive[i].end(), CBigNum(0));
        sPolyNegative[i] = circuit.wCj[i-1-m];
    }

    rPolyPositive[2*m+1] = D;
    fill(rPolyNegative[2*m+1].begin(), rPolyNegative[2*m+1].end(), CBigNum(0));
    fill(sPolyNegative[2*m+1].begin(), sPolyNegative[2*m+1].end(), CBigNum(0));

    // rDashPoly
    CBN_matrix rDashPolyPositive(2*m+2, CBN_vector(n));
    CBN_matrix rDashPolyNegative(2*m+2, CBN_vector(n));

    fill(rDashPolyPositive[0].begin(), rDashPolyPositive[0].end(), CBigNum(0));
    fill(rDashPolyNegative[0].begin(), rDashPolyNegative[0].end(), CBigNum(0));


    for(int i=1; i<2*m+1; i++)
        for(int j=0; j<n; j++) {
            rDashPolyPositive[i][j] = rPolyPositive[i][j].mul_mod(circuit.YDash[j],q);
            rDashPolyPositive[i][j] = (rDashPolyPositive[i][j] + 2 * sPolyPositive[i][j]) % q;
            rDashPolyNegative[i][j] = rPolyNegative[i][j].mul_mod(circuit.YDash[j],q);
            rDashPolyNegative[i][j] = (rDashPolyNegative[i][j] + 2 * sPolyNegative[i][j]) % q;
        }

    for(int j=0; j<n; j++)
            rDashPolyPositive[2*m+1][j] = D[j].mul_mod(circuit.YDash[j],q);

    fill(rDashPolyNegative[2*m+1].begin(), rDashPolyNegative[2*m+1].end(), CBigNum(0));

    // tPoly
    CBN_vector tPoly(7*m+3);
    CBigNum tcoef;
    CBN_vector *oper1, *oper2;

    for(int k=0; k<7*m+3; k++) {
        tcoef = CBigNum(0);
        for(int i=max(k-5*m-1,-m); i<min(k-m,2*m+1)+1; i++) {
            int j = k - 3*m - i;
            oper1 = i > 0 ? &rPolyPositive[i] : &rPolyNegative[-i];
            oper2 = j > 0 ? &rDashPolyPositive[j] : &rDashPolyNegative[-j];
            tcoef += dotProduct(*oper1, *oper2, q);
            tcoef %=q;
        }
        tPoly[k] = tcoef;
    }

    // sanity check
    if (tPoly[3*m] != (2*circuit.Kconst)%q)
        throw std::runtime_error("SerialNumberSoK_small - error: sanity check failed");

    tPoly[3*m] = CBigNum(0);

    // commit to the polynomial
    polyComm.Commit(tPoly);


    // ****************************************************************************
    // *********************** STEP 4: Challenge Component ************************
    // ****************************************************************************

    hasher << polyComm.U.ToString();
    for(unsigned int i=0; i<m1dash; i++) hasher << polyComm.Tf[i].ToString();
    for(unsigned int i=0; i<m1dash; i++) hasher << polyComm.Trho[i].ToString();

    // get the challenge component x
    CBigNum x = CBigNum(hasher.GetHash()) % q;

    // precomputation of x powers
    CBN_vector xPowersPos(m2dash*ndash+1);
    CBN_vector xPowersNeg(m1dash*ndash+1);
    xPowersPos[0] = xPowersNeg[0] = CBigNum(1);
    xPowersPos[1] = x;
    xPowersNeg[1] = x.pow_mod(-1,q);
    for(int i=2; i<m2dash*ndash+1; i++)
        xPowersPos[i] = xPowersPos[i-1].mul_mod(x,q);
    for(int i=2; i<m1dash*ndash+1; i++)
        xPowersNeg[i] = xPowersNeg[i-1].mul_mod(xPowersNeg[1],q);


    // ****************************************************************************
    // **************************** STEP 5: Poly Eval *****************************
    // ****************************************************************************

    // evaluate the polynomial at x
    polyComm.Eval(xPowersPos, xPowersNeg);

    CBN_vector r_vec(n+pads, CBigNum(0));

    for(unsigned int j=0; j<n; j++){
        for(unsigned int rcoef=0; rcoef<2*m+2; rcoef++) {
            r_vec[j] +=
                    rPolyPositive[rcoef][j].mul_mod(xPowersPos[rcoef],q) +
                    rPolyNegative[rcoef][j].mul_mod(xPowersNeg[rcoef],q);
            r_vec[j] %= q;
        }
    }

    CBN_vector s_vec(n+pads, CBigNum(0));

    for(unsigned int j=0; j<n; j++){
        for(unsigned int scoef=0; scoef<2*m+2; scoef++) {
            s_vec[j] +=
                    sPolyPositive[scoef][j].mul_mod(xPowersPos[scoef],q) +
                    sPolyNegative[scoef][j].mul_mod(xPowersNeg[scoef],q);
            s_vec[j] %= q;
        }
    }

    rho = f_delta.mul_mod(xPowersPos[2*m+1],q);

    for(unsigned int i=1; i<m+1; i++) {
        rho +=
                f_alpha[i-1].mul_mod(xPowersPos[i].mul_mod(y.pow_mod(i,q),q),q) +
                f_beta[i-1].mul_mod(xPowersNeg[i],q) +
                f_gamma[i-1].mul_mod(xPowersPos[m+i],q);
        rho %= q;
    }

    // ****************************************************************************
    // ********************** STEP 6: Inner Product Argument **********************
    // ****************************************************************************

    CBN_vector y_vec;
    CBigNum ym = y.pow_mod(m, q);
    CBigNum temp = CBigNum(1);
    for(unsigned int j=0; j<n+pads; j++) {
        temp = temp.mul_mod(ym,q);
        y_vec.push_back(temp);
    }

    CBN_vector temp_vec = hadamard(y_vec, r_vec, q);
    CBN_vector temp_vec2;
    for(unsigned j=0; j<s_vec.size(); j++)
        temp_vec2.push_back(s_vec[j].mul_mod(CBigNum(2), q));

    CBN_vector rdash_vec1(n+pads);
    for(unsigned int j=0; j<n+pads; j++)
        rdash_vec1[j] = (temp_vec[j] + temp_vec2[j]) % q;

    y_vec.clear();
    ym = y.pow_mod(-m, q);
    temp = CBigNum(1);
    for(unsigned int j=0; j<n+pads; j++) {
        temp = temp.mul_mod(ym,q);
        y_vec.push_back(temp.mul_mod(CBigNum(2),q));
    }
    temp_vec.clear();
    temp_vec = hadamard(y_vec, s_vec, q);

    CBN_vector rdash_vec2(n+pads);
    for(unsigned int j=0; j<n+pads; j++)
        rdash_vec2[j] = (r_vec[j] + temp_vec[j]) % q;


    // Inner-product PROVE
    CBigNum ComR = pedersenCommitment(params, r_vec, CBigNum(0));
    comRdash = pedersenCommitment(params, rdash_vec2, CBigNum(0));

    CBigNum Pinner = ComR.mul_mod(comRdash, p);

    CBigNum z = dotProduct(r_vec, rdash_vec1, q);

    CBN_matrix ck_inner_g = ck_inner_gen(params);

    CBN_matrix r1(1, CBN_vector(r_vec));
    CBN_matrix r2(1, CBN_vector(rdash_vec1));
    innerProduct.Prove(ck_inner_g, Pinner, z, r1, r2, y);

    // Remove y1 from ComC
    ComC.pop_back();
}


bool SerialNumberSoK_small::Verify(const CBigNum& coinSerialNumber,
        const CBigNum& valueOfCommitmentToCoin, const uint256 msghash) const
{
    std::vector<SerialNumberSoKProof> proof(1, SerialNumberSoKProof(*this, coinSerialNumber, valueOfCommitmentToCoin, msghash));

    return SerialNumberSoKProof::BatchVerify(proof);

}


bool SerialNumberSoKProof::BatchVerify(std::vector<SerialNumberSoKProof> &proofs) {
    const CBigNum q = proofs[0].signature.params->serialNumberSoKCommitmentGroup.groupOrder;
    const CBigNum p = proofs[0].signature.params->serialNumberSoKCommitmentGroup.modulus;
    const int m =  ZKP_M;
    const int n =  ZKP_N;
    const int m1dash =  ZKP_M1DASH;
    const int m2dash =  ZKP_M2DASH;
    const int ndash =  ZKP_NDASH;
    const int pads = ZKP_PADS;

    // ****************************************************************************
    // **************************** STEP 1: Parsing *******************************
    // ****************************************************************************

    uint256 msghash;
    CBigNum S;
    CBigNum y1;
    CBN_vector ComA;
    CBN_vector ComB;
    CBN_vector ComC;
    CBigNum ComD;
    CBigNum comRdash;
    CBigNum rho;
    PolynomialCommitment *polyComm;
    Bulletproofs *innerProduct;
    const ZerocoinParams *params;

    CBigNum ny;
    CBN_vector ymPowers;
    CBigNum temp;

    for(unsigned int i=0; i<proofs.size(); i++)
    {
        msghash = proofs[i].msghash;
        S = proofs[i].coinSerialNumber;
        y1 = proofs[i].valueOfCommitmentToCoin;
        ComA = proofs[i].signature.ComA;
        ComB = proofs[i].signature.ComB;
        ComC = proofs[i].signature.ComC;
        ComD = proofs[i].signature.ComD;
        comRdash  = proofs[i].signature.comRdash;
        rho = proofs[i].signature.rho;
        polyComm = &proofs[i].signature.polyComm;
        innerProduct = &proofs[i].signature.innerProduct;

        // Restore y1 in ComC
        CBN_vector ComC_(ComC);
        ComC_.push_back(y1);

        // Assert inputs in correct groups
        if( S < CBigNum(0) || S > CBigNum(2).pow(256))
            throw std::runtime_error("wrong value for S");

        if( ComD < CBigNum(0) || ComD > p )
            throw std::runtime_error("wrong value for ComD");

        for(int i=0; i<m; i++) {
            if( ComA[i] < CBigNum(0) || ComA[i] > p )
                throw std::runtime_error("wrong value for ComA at " + to_string(i));
            if( ComB[i] < CBigNum(0) || ComB[i] > p )
                throw std::runtime_error("wrong value for ComB at " + to_string(i));
            if( ComC_[i] < CBigNum(0) || ComC_[i] > p )
                throw std::runtime_error("wrong value for ComC at " + to_string(i));
        }

        if( comRdash < CBigNum(0) || comRdash > p )
            throw std::runtime_error("wrong value for comRdash");

        for(int i=0; i<m1dash; i++)
            if( polyComm->Tf[i] < CBigNum(0) || polyComm->Tf[i] > p )
                throw std::runtime_error("wrong value for Tf at " + to_string(i));

        for(int i=0; i<m2dash; i++)
            if( polyComm->Trho[i] < CBigNum(0) || polyComm->Trho[i] > p )
                throw std::runtime_error("wrong value for Trho at " + to_string(i));

        if( polyComm->U < CBigNum(0) || polyComm->U > p )
            throw std::runtime_error("wrong value for U");

        for(int i=0; i<ndash; i++)
            if( polyComm->tbar[i] < CBigNum(0) || polyComm->tbar[i] > q )
                throw std::runtime_error("wrong value for tbar at " + to_string(i));

        if( polyComm->taubar < CBigNum(0) || polyComm->taubar > q )
            throw std::runtime_error("wrong value for taubar");

        for(int j=0; j<(int)innerProduct->pi[0].size(); j++)
            if( innerProduct->pi[0][j] < CBigNum(0) || innerProduct->pi[0][j] > p )
                throw std::runtime_error("wrong value for pi[0] at j=" + to_string(j));

        for(int j=0; j<(int)innerProduct->pi[1].size(); j++)
            if( innerProduct->pi[1][j] < CBigNum(0) || innerProduct->pi[1][j] > p )
                throw std::runtime_error("wrong value for pi[1] at j=" + to_string(j));

        const int M1 = innerProduct->final_a.size();
        const int N1 = innerProduct->final_a[0].size();

        for(int i=0; i<M1; i++)
            for(int j=0; j<N1; j++) {
                if( innerProduct->final_a[i][j] < CBigNum(0) || innerProduct->final_a[i][j] > q )
                    throw std::runtime_error("wrong value for final_a at " + to_string(i) + " " + to_string(j));
                if( innerProduct->final_b[i][j] < CBigNum(0) || innerProduct->final_b[i][j] > q )
                    throw std::runtime_error("wrong value for final_b at " + to_string(i) + " " + to_string(j));
            }


    }


    // ****************************************************************************
    // *********************** STEP 2: Compute Challenges *************************
    // ****************************************************************************

    std::vector<DelayedProof> proofsForLater;

    for(unsigned int i=0; i<proofs.size(); i++)
    {
        msghash = proofs[i].msghash;
        S = proofs[i].coinSerialNumber;
        y1 = proofs[i].valueOfCommitmentToCoin;
        ComA = proofs[i].signature.ComA;
        ComB = proofs[i].signature.ComB;
        ComC = proofs[i].signature.ComC;
        ComD = proofs[i].signature.ComD;
        comRdash  = proofs[i].signature.comRdash;
        rho = proofs[i].signature.rho;
        polyComm = &proofs[i].signature.polyComm;
        innerProduct = &proofs[i].signature.innerProduct;
        params = proofs[i].signature.params;

        // Restore y1 in ComC
        CBN_vector ComC_(ComC);
        ComC_.push_back(y1);

        CHashWriter1024 hasher(0,0);
        hasher << msghash << ComD.ToString();

        for(int i=0; i<m; i++)
            hasher << ComA[i].ToString() << ComB[i].ToString() << ComC_[i].ToString();

        // get the challenge component y
        CBigNum y = CBigNum(hasher.GetHash()) % q;

        ny = y.pow_mod(-ZKP_M, q);
        ymPowers.clear();
        ymPowers.push_back(CBigNum(1));
        temp = CBigNum(1);

        for(unsigned int i=0; i<n+pads; i++) {
            temp = temp.mul_mod(ny, q);
            ymPowers.push_back(temp);
        }

        hasher << polyComm->U.ToString();
        for(int i=0; i<m1dash; i++) hasher << polyComm->Tf[i].ToString();
        for(int i=0; i<m1dash; i++) hasher << polyComm->Trho[i].ToString();

        // get the challenge component x
        CBigNum x = CBigNum(hasher.GetHash()) % q;

        // precomputation of x powers
        CBN_vector xPowersPos(m2dash*ndash+1);
        CBN_vector xPowersNeg(m1dash*ndash+1);
        xPowersPos[0] = xPowersNeg[0] = CBigNum(1);
        xPowersPos[1] = x;
        xPowersNeg[1] = x.pow_mod(-1,q);
        for(int i=2; i<m2dash*ndash+1; i++)
            xPowersPos[i] = xPowersPos[i-1].mul_mod(x,q);
        for(int i=2; i<m1dash*ndash+1; i++)
            xPowersNeg[i] = xPowersNeg[i-1].mul_mod(xPowersNeg[1],q);



        // ****************************************************************************
        // ********************** STEP 3: Find rDash and Kconst ***********************
        // ****************************************************************************

        // set arithmetic circuit wire constraints
        ArithmeticCircuit circuit(params);
        circuit.setConstraints(S);

        // set circuit w-Polynomials
        circuit.setYPoly(y);

        // s_vec
        CBN_vector s_vec(n);
        CBigNum exp;
        CBigNum yneg = y.pow_mod(CBigNum(-1),q);

        for(int j=0; j<n; j++){
            exp = CBigNum(1);
            s_vec[j] = CBigNum(0);
            for(int i=1; i<m+1; i++) {
                exp = exp.mul_mod(yneg, q);
                s_vec[j] +=
                        circuit.wAj[i-1][j].mul_mod(exp.mul_mod(xPowersNeg[i],q),q) +
                        circuit.wBj[i-1][j].mul_mod(xPowersPos[i],q) +
                        circuit.wCj[i-1][j].mul_mod(xPowersNeg[m+i],q);
                s_vec[j] %= q;
            }
            s_vec[j] *= 2;
            s_vec[j] %= q;
        }

        // ****************************************************************************
        // ******************************* FINAL STEP *********************************
        // ****************************************************************************

        // ComR
        CBN_vector zerovec(1, CBigNum(0));
        CBigNum ComR = pedersenCommitment(params, zerovec, -rho);

        ComR = ComR.mul_mod(ComD.pow_mod(xPowersPos[2*m+1],p),p);

        for(int i=1; i<m+1; i++) {
            ComR = ComR.mul_mod(ComA[i-1].pow_mod(xPowersPos[i].mul_mod(y.pow_mod(i,q),q),p),p);
            ComR = ComR.mul_mod(ComB[i-1].pow_mod(xPowersNeg[i],p),p);
            ComR = ComR.mul_mod(ComC_[i-1].pow_mod(xPowersPos[m+i],p),p);
        }

        // restore PolynomialCommitment object from commitments
        PolynomialCommitment polyCommitment(params, polyComm->Tf, polyComm->Trho, polyComm->U,
                polyComm->tbar, polyComm->taubar, xPowersPos, xPowersNeg);

        // verify the polynomial commitment and save the evaluation in z
        CBigNum z;
        if(!polyCommitment.Verify(z))
            return false;

        // verify the inner product argument
        z = (z + 2*circuit.Kconst) % q;

        DelayedProof dp(innerProduct);
        dp.ComR = ComR;
        dp.comRdash = comRdash;
        dp.z = z;
        dp.ymPowers = ymPowers;
        dp.s_vec2 = s_vec;

        proofsForLater.push_back(dp);
    }

    params = proofs[0].signature.params;

    unsigned int vec2size = proofsForLater[0].s_vec2.size();
    CBN_vector test_vec(vec2size, CBigNum(0));
    CBigNum comTest = CBigNum(1);

    std::vector<DelayedProof> inner_product_proof;
    CBigNum gamma;
    DelayedProof *proof;

    for(unsigned int k=0; k<proofsForLater.size(); k++) {
        proof = &proofsForLater[k];
        inner_product_proof.push_back(*proof);

        gamma = CBigNum::randBignum(q);

        for(unsigned int i=0; i<vec2size; i++)
            test_vec[i] = (test_vec[i] + gamma.mul_mod(proof->s_vec2[i].mul_mod(proof->ymPowers[i+1],q),q)) % q;

        comTest = comTest.mul_mod(
                ((proof->ComR.pow_mod(-1,p)).mul_mod(proof->comRdash,p)).pow_mod(gamma,p),p);
    }

    CBigNum test = pedersenCommitment(params, test_vec, CBigNum(0));

    if(test != comTest) {
        cout << "\ndifferent test and comTest" << endl;
        cout << "test = " << test.ToString() << endl;
        cout << "comTest = " << comTest.ToString() << endl;
        return false;
    }

    CBN_matrix ck_inner_g = ck_inner_gen(params);
    bool valid = BatchBulletproofs(params, ck_inner_g, inner_product_proof);

    return valid;

}

bool SerialNumberSoKProof::BatchBulletproofs(const ZerocoinParams* ZCp, const CBN_matrix ck_inner_g, std::vector<DelayedProof> &proofs) {
    const CBigNum q = ZCp->serialNumberSoKCommitmentGroup.groupOrder;
    const CBigNum p = ZCp->serialNumberSoKCommitmentGroup.modulus;
    const CBigNum u_inner_prod = ZCp->serialNumberSoKCommitmentGroup.u_inner_prod;

    // Initialize
    DelayedProof dp = proofs[0];
    int N1 = dp.innerProduct.pi[0].size();

    CBigNum Ptest = CBigNum(1);

    std::vector<fBE> forBigExpo;
    CBigNum gamma, x1, u_inner, P_inner;
    CBigNum x, Ak, Bk;
    CBN_vector xlist;
    CBigNum pt1, pt2;
    CBigNum A, B, z;
    for(unsigned int i=0; i<proofs.size(); i++)
    {
        dp = proofs[i];
        A = dp.ComR;
        B = dp.comRdash;
        z = dp.z;

        gamma = CBigNum::randBignum(q);

        CBigNum P_inner_prod = A.mul_mod(B, p);

        // Inserting the z into u
        CHashWriter1024 hasher(0,0);
        hasher << u_inner_prod.ToString() << P_inner_prod.ToString() << z.ToString();
        x1 = CBigNum(hasher.GetHash()) % q;

        u_inner = u_inner_prod.pow_mod(x1,p);
        P_inner = P_inner_prod.mul_mod(u_inner.pow_mod(z,p),p);

        // Starting the actual protocol
        xlist.clear();

        for(int i=0; i<N1; i++) {
            Ak = dp.innerProduct.pi[0][i];
            Bk = dp.innerProduct.pi[1][i];

            hasher << Ak.ToString() << Bk.ToString();
            x = CBigNum(hasher.GetHash()) % q;

            xlist.push_back(x);

            P_inner = P_inner.mul_mod(
                    (Ak.pow_mod(x.pow_mod(2,q),p)).mul_mod(Bk.pow_mod(x.pow_mod(-2,q),p),p),p);
        }

        z = dp.innerProduct.final_a[0][0].mul_mod(dp.innerProduct.final_b[0][0],q);

        pt1 = P_inner.pow_mod(gamma,p);
        pt2 = u_inner.pow_mod(z.mul_mod(-gamma,q),p);
        Ptest = Ptest.mul_mod( pt1.mul_mod(pt2,p) ,p);

        fBE new_element;
        new_element.gamma = gamma;
        new_element.xlist = xlist;
        new_element.ymPowers = dp.ymPowers;
        new_element.a = dp.innerProduct.final_a[0][0];
        new_element.b = dp.innerProduct.final_b[0][0];

        forBigExpo.push_back(new_element);
    }

    CBN_vector gh_final = getFinal_gh(ZCp, ck_inner_g[0], forBigExpo);

    return (gh_final[0].mul_mod(gh_final[1],p) == Ptest);
}


CBN_vector SerialNumberSoKProof::getFinal_gh(const ZerocoinParams* ZCp, CBN_vector gs, std::vector<fBE> forBigExpo)
{
    const CBigNum q = ZCp->serialNumberSoKCommitmentGroup.groupOrder;
    const CBigNum p = ZCp->serialNumberSoKCommitmentGroup.modulus;

    int logn = forBigExpo[0].xlist.size();
    int n = gs.size();
    CBN_vector sg_expo(n, CBigNum(0));
    CBN_vector sh_expo(n, CBigNum(0));
    CBN_vector xnlist;

    for(int k=0; k<(int)forBigExpo.size(); k++) {
        fBE comp = forBigExpo[k];

        std::reverse(comp.xlist.begin(),comp.xlist.end());

        xnlist.clear();
        for(int i=0; i<logn; i++)
            xnlist.push_back(comp.xlist[i].pow_mod(-1,q));

        std::vector< std::vector<int>> binary_lookup = Bulletproofs::findBinaryLookup(logn);

        CBN_vector temp_g, temp_h;
        CBigNum sg_i, sh_i;
        std::vector<int> bi;
        for(int i=0; i<n; i++) {
            sg_i = (comp.gamma).mul_mod(comp.a,q);
            sh_i = (comp.gamma).mul_mod((comp.b).mul_mod(comp.ymPowers[i+1],q),q);
            bi = binary_lookup[i];

            for(int j=0; j<logn; j++) {
                if (bi[j] == 1) {
                    sg_i = sg_i.mul_mod(comp.xlist[j],q);
                    sh_i = sh_i.mul_mod(xnlist[j],q);
                } else {
                    sg_i = sg_i.mul_mod(xnlist[j],q);
                    sh_i = sh_i.mul_mod(comp.xlist[j],q);
                }
            }

            temp_g.push_back(sg_i);
            temp_h.push_back(sh_i);

            sg_expo[i] = (sg_expo[i] + sg_i) % q;
            sh_expo[i] = (sh_expo[i] + sh_i) % q;
        }
    }

    CBN_vector gh_final(2, CBigNum(1));

    for(int i=0; i<n; i++) {
        gh_final[0] = gh_final[0].mul_mod(gs[i].pow_mod(sg_expo[i],p),p);
        gh_final[1] = gh_final[1].mul_mod(gs[i].pow_mod(sh_expo[i],p),p);
    }

    return gh_final;
}

} /* namespace libzerocoin */


