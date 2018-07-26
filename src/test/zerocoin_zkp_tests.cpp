#include "chainparams.h"
#include "libzerocoin/ArithmeticCircuit.h"
#include "libzerocoin/PolynomialCommitment.h"
#include "libzerocoin/InnerProductArgument.h"
#include "libzerocoin/SerialNumberSoK_small.h"
#include "libzerocoin/SerialNumberSignatureOfKnowledge.h"
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <time.h>


using namespace libzerocoin;

BOOST_AUTO_TEST_SUITE(zerocoin_zkp_tests)

/*
BOOST_AUTO_TEST_CASE(parameters_tests)
{
    std::cout << endl;
    std::cout << "*** parameters_tests ***" << endl;
    std::cout << "------------------------" << endl;

    SelectParams(CBaseChainParams::MAIN);
    ZerocoinParams *ZCParams = Params().Zerocoin_Params(false);
    (void)ZCParams;

    std::cout << "- Testing generators..." << endl;
    for(unsigned int i=0; i<512; i++)
        BOOST_CHECK_MESSAGE( ZCParams->serialNumberSoKCommitmentGroup.gis[i].pow_mod(
                ZCParams->serialNumberSoKCommitmentGroup.groupOrder,
                ZCParams->serialNumberSoKCommitmentGroup.modulus) == CBigNum(1),
                "Generator gis[i] error for i=" << i << "\n");

    std::cout << endl;
}
*/

BOOST_AUTO_TEST_CASE(arithmetic_circuit_tests)
{
    std::cout << "*** arithmetic_circuit_tests ***" << endl;
    std::cout << "--------------------------------" << endl;

    SelectParams(CBaseChainParams::MAIN);
    ZerocoinParams *ZCParams = Params().Zerocoin_Params(false);
    (void)ZCParams;

    CBigNum a = ZCParams->coinCommitmentGroup.g;
    CBigNum b = ZCParams->coinCommitmentGroup.h;
    CBigNum q = ZCParams->serialNumberSoKCommitmentGroup.groupOrder;

    // mint a coin
    PrivateCoin coin(ZCParams, CoinDenomination::ZQ_ONE);
    // get random Y
    CBigNum Y = CBigNum::randBignum(q);

    ArithmeticCircuit circuit(ZCParams);
    circuit.setWireValues(coin);
    circuit.setYPoly(Y);

    // If multiplication gates hold this should be true
    std::cout << "- Testing A times B equals C..." << endl;
    for(unsigned int i=0; i<ZKP_M; i++) for(unsigned int j=0; j<ZKP_N; j++) {
        BOOST_CHECK_MESSAGE(
                circuit.A[i][j].mul_mod(circuit.B[i][j], q) == circuit.C[i][j],
                "Circuit Specification:: Hadamard Test failed\n" <<
                "(i=" << i << ", j=" << j << ")\n" <<
                "A[i][j]=" << circuit.A[i][j] << "\n" <<
                "B[i][j]=" << circuit.B[i][j] << "\n" <<
                "C[i][j]=" << circuit.C[i][j] << "\n");
    }

    // If circuit correctly evaluates (a^serial)*(b^randomness) this should be true
    std::cout << "- Testing C_final equals Logarithm..." << endl;
    CBigNum logarithm =
            a.pow_mod(circuit.getSerialNumber(),q).mul_mod(
            b.pow_mod(circuit.getRandomness(),q),q);
    CBigNum Cfinal = circuit.C[ZKP_M-1][0];
    BOOST_CHECK_MESSAGE( logarithm == Cfinal,
            "Circuit Specification:: Correctness Test failed\n" <<
            "logarithm = " << logarithm.ToString() << "\n" <<
            "Cfinal = " << Cfinal.ToString() << "\n");


    // Checking that the expressions in Equation (1) of the paper hold
    std::cout << "- Testing the Arithmetic Constraints (eq. 1)..." << endl;
    for(unsigned int i=0; i<4*ZKP_SERIALSIZE-2; i++) {
        BOOST_CHECK_MESSAGE( circuit.sumWiresDotWs(i) == (circuit.K[i] % q),
                "Circuit Specification:: Arithmetic Constraints Test failed at i=" << i << "\n" <<
                "sumWiresDotWs(i) = " << circuit.sumWiresDotWs(i) << "\n" <<
                "K[i] = " << circuit.K[i].ToString() << "\n");
    }

    // Checking that the expressions in Equation (2) of the paper hold
    //std::cout << "- Testing Polynomials Evaluate (eq. 2)..." << endl;
    //CBigNum sum = circuit.sumWiresDotWPoly();
    //BOOST_CHECK_MESSAGE( sum == circuit.Kconst,
    //        "Circuit Specification:: Constraints Polynomial Test failed\n" <<
    //        "sum = " << sum.ToString() << "\n" <<
    //        "Kconst = " << circuit.Kconst << "\n");

    // New circuit with random assignment
    ArithmeticCircuit newCircuit(circuit);
    for(unsigned int i=0; i<ZKP_M; i++) {
        random_vector_mod(newCircuit.A[i], q);
        random_vector_mod(newCircuit.B[i], q);
        for(unsigned int j=0; j<ZKP_N; j++)
            newCircuit.C[i][j] = newCircuit.A[i][j].mul_mod(newCircuit.B[i][j],q);
    }

    // If circuit correctly evaluates (a^serial)*(b^randomness) we have a problem
    std::cout << "- Testing C_final != Logarithm for wrong assignment..." << endl;
    logarithm =
            a.pow_mod(newCircuit.getSerialNumber(),q).mul_mod(
            b.pow_mod(newCircuit.getRandomness(),q),q);
    Cfinal = newCircuit.C[ZKP_M-1][0];
    BOOST_CHECK_MESSAGE( logarithm != Cfinal,
            "Circuit Specification:: Correctness Test passed for wrong assignment\n");

    // Checking that the expressions in Equation (1) of the paper does not hold
    std::cout << "- Testing the Arithmetic Constraints (eq. 1) for wrong assignment..." << endl;
    for(unsigned int i=0; i<4*ZKP_SERIALSIZE-2; i++) {
        BOOST_CHECK_MESSAGE( newCircuit.sumWiresDotWs(i) != newCircuit.K[i],
                "Circuit Specification:: Arithmetic Constraints Test passed at i=" << i << " with wrong assignment");
    }

    // Checking that the expressions in Equation (2) of the does not paper hold
    //std::cout << "- Testing Polynomials Evaluate (eq. 2) for wrong assignment..." << endl;
    //BOOST_CHECK_MESSAGE( newCircuit.sumWiresDotWPoly() != newCircuit.Kconst,
    //        "Circuit Specification:: Constraints Polynomial Test passed with wrong assignment\n");

    std::cout << endl;
}

/*
// Evaluate tpolynomial at x
CBigNum eval_tpoly(CBN_vector tpoly, CBN_vector xPowersPos, CBN_vector xPowersNeg, CBigNum q)
{
    CBigNum sum = CBigNum(0);
    for(unsigned int i=0; i<=ZKP_NDASH*ZKP_M1DASH; i++)
        sum = ( sum + tpoly[i].mul_mod(xPowersNeg[ZKP_NDASH*ZKP_M1DASH-i],q) ) % q;
    for(unsigned int i=ZKP_NDASH*ZKP_M1DASH+1; i<=ZKP_NDASH*(ZKP_M1DASH+ZKP_M2DASH); i++)
        sum = ( sum + tpoly[i].mul_mod(xPowersPos[i-ZKP_NDASH*ZKP_M1DASH],q) ) % q;
    return sum;
}

BOOST_AUTO_TEST_CASE(polynomial_commitment_tests)
{
    std::cout << "*** polynomial_commitment_tests ***" << endl;
    std::cout << "-----------------------------------" << endl;

    SelectParams(CBaseChainParams::MAIN);
    ZerocoinParams *ZCParams = Params().Zerocoin_Params(false);
    (void)ZCParams;

    CBigNum q = ZCParams->serialNumberSoKCommitmentGroup.groupOrder;
    CBigNum p = ZCParams->serialNumberSoKCommitmentGroup.modulus;

    // generate a random tpolynomial with 0 constant term
    CBN_vector tpoly(ZKP_NDASH*(ZKP_M1DASH+ZKP_M2DASH)+1);
    for(unsigned int i=0; i<tpoly.size(); i++)
        tpoly[i] = CBigNum::randBignum(q);
    tpoly[ZKP_M1DASH*ZKP_NDASH] = CBigNum(0);

    // generate a random evaluation point x in R and compute powers
    CBigNum x = CBigNum::randBignum(q);
    CBN_vector xPowersPositive(ZKP_M2DASH*ZKP_NDASH+1);
    CBN_vector xPowersNegative(ZKP_M1DASH*ZKP_NDASH+1);
    xPowersPositive[0] = xPowersNegative[0] = CBigNum(1);
    xPowersPositive[1] = x;
    xPowersNegative[1] = x.pow_mod(-1,q);
    for(unsigned int i=2; i<ZKP_M2DASH*ZKP_NDASH+1; i++)
        xPowersPositive[i] = x.pow_mod(i,q);
    for(unsigned int i=2; i<ZKP_M1DASH*ZKP_NDASH+1; i++)
        xPowersNegative[i] = x.pow_mod(-(int)i,q);

    // Poly-Commit and Poly-Evaluate
    PolynomialCommitment polyCommitment(ZCParams);
    polyCommitment.Commit(tpoly);
    polyCommitment.Eval(xPowersPositive, xPowersNegative);


    // Poly-Verify: For honest prover, verifier should be satisfied
    std::cout << "- Testing PolyVerify for honest prover..." << endl;
    CBigNum val;    // val = t(x) if proofs checks out
    BOOST_CHECK_MESSAGE( polyCommitment.Verify(val),
            "Polynomial Commitment:: PoliVerify returned FALSE.\n");

    // Poly-Verify: For honest prover, verifier is able to compute t(x)
    std::cout << "- Testing t(x) == dotProduct(tbar,xPowersPos)..." << endl;
    CBigNum tx = eval_tpoly(tpoly, xPowersPositive, xPowersNegative, q);
    BOOST_CHECK_MESSAGE( val == tx,
            "Polynomial Commitment:: Verifier computed wrong t(x) value.\n" <<
            "t(x) = " << tx.ToString() << "\n" <<
            "<tbar, [1,...,x^n]> = " << val.ToString() << "\n");

    // Create copies of the polynomial commitment and mess things up
    PolynomialCommitment newPolyComm1(polyCommitment);
    PolynomialCommitment newPolyComm2(polyCommitment);
    PolynomialCommitment newPolyComm3(polyCommitment);
    random_vector_mod(newPolyComm1.tbar, q);
    random_vector_mod(newPolyComm2.Tf, q);
    random_vector_mod(newPolyComm3.Trho, q);

    // Poly-Verify: For dishonest prover, verifier should fail the test
    std::cout << "- Testing PolyVerify for dishonest prover..." << endl;
    BOOST_CHECK_MESSAGE( !newPolyComm1.Verify(val),
            "Polynomial Commitment:: PoliVerify returned TRUE[1] for dishonest prover.\n");
    BOOST_CHECK_MESSAGE( !newPolyComm2.Verify(val),
            "Polynomial Commitment:: PoliVerify returned TRUE[2] for dishonest prover.\n");
    BOOST_CHECK_MESSAGE( !newPolyComm3.Verify(val),
            "Polynomial Commitment:: PoliVerify returned TRUE[3] for dishonest prover.\n");

    std::cout << endl;
}

BOOST_AUTO_TEST_CASE(inner_product_argument_tests)
{
    std::cout << "*** inner_product_argument_tests ***" << endl;
    std::cout << "------------------------------------" << endl;

    SelectParams(CBaseChainParams::MAIN);
    ZerocoinParams *ZCParams = Params().Zerocoin_Params(false);
    (void)ZCParams;

    CBigNum q = ZCParams->serialNumberSoKCommitmentGroup.groupOrder;
    CBigNum p = ZCParams->serialNumberSoKCommitmentGroup.modulus;

    // Get random y in Z_q and (N+PADS)-vectors
    CBigNum y = CBigNum::randBignum(q);
    CBN_vector a_sets(ZKP_N+ZKP_PADS);
    CBN_vector b_sets(ZKP_N+ZKP_PADS);
    random_vector_mod(a_sets, q);
    random_vector_mod(b_sets, q);


    // Inner-product PROVE
    InnerProductArgument innerProduct(ZCParams);
    innerProduct.Prove(y, a_sets, b_sets);

    // compute ck_inner sets
    pair<CBN_vector, CBN_vector> resultSets = innerProduct.ck_inner_gen(ZCParams, y);
    CBN_vector ck_inner_g = resultSets.first;
    CBN_vector ck_inner_h = resultSets.second;

    // Compute commitment A to a_sets under ck_inner_g
    CBigNum A = CBigNum(1);
    for(unsigned int i=0; i<a_sets.size(); i++)
        A = A.mul_mod(ck_inner_g[i].pow_mod(a_sets[i], p), p);

    // Compute commitment B to b_sets under ck_inner_h
    CBigNum B = CBigNum(1);
    for(unsigned int i=0; i<b_sets.size(); i++)
        B = B.mul_mod(ck_inner_h[i].pow_mod(b_sets[i], p), p);

    // Inner product z = <a,b>
    CBigNum z = dotProduct(a_sets, b_sets, q);


    // random_z != z
    CBigNum random_z = CBigNum::randBignum(q);
    while( random_z == dotProduct(a_sets, b_sets, q) ) random_z = CBigNum::randBignum(q);

    // innerProductVerify
    std::cout << "- Testing innerProductVerify..." << endl;
    bool res = innerProduct.Verify(ZCParams, y, A, B, z);

    BOOST_CHECK_MESSAGE(res,"InnerProduct:: Verification failed\n");

    // Inner product z != <a, b>  (z = 0, z = 1, z = random)
    std::cout << "- Testing innerProductVerify for dishonest prover..." << endl;

    BOOST_CHECK_MESSAGE( !innerProduct.Verify(ZCParams, y, A, B, CBigNum(0)),
            "InnerProduct:: Verification returned TRUE[1] for dishonest prover\n");
    BOOST_CHECK_MESSAGE( !innerProduct.Verify(ZCParams, y, A, B, CBigNum(1)),
            "InnerProduct:: Verification returned TRUE[2] for dishonest prover\n");
    BOOST_CHECK_MESSAGE( !innerProduct.Verify(ZCParams, y, A, B, random_z),
            "InnerProduct:: Verification returned TRUE[3] for dishonest prover\n");


    std::cout << endl;
}
*/

BOOST_AUTO_TEST_CASE(signature_of_knowledge_tests)
{
    std::cout << "*** signature_of_knowledge_tests DEBUG ***" << endl;
    std::cout << "------------------------------------------" << endl;

    SelectParams(CBaseChainParams::MAIN);
    ZerocoinParams *ZCParams = Params().Zerocoin_Params(false);
    (void)ZCParams;

    // create message hash
    CHashWriter hasher(0,0);
    hasher << std::string("Chancellor on brink of second bailout for banks");
    uint256 msghash = hasher.GetHash();

    // mint a coin
    PrivateCoin newCoin(ZCParams, CoinDenomination::ZQ_ONE);

    // commit to this coin
    const CBigNum newCoin_value = newCoin.getPublicCoin().getValue();
    Commitment commitment(&(ZCParams->serialNumberSoKCommitmentGroup), newCoin_value);


    std::cout << "- Creating the Signature of Knowledge..." << endl;

    // create the signature of knowledge
    SerialNumberSoK_small sigOfKnowledge(ZCParams, newCoin, commitment, msghash);

    std::cout << "- Serializing the Signature of Knowledge..." << endl;

    // serialize the SoK to a CDataStream object.
    CDataStream serializedSoK(SER_NETWORK, PROTOCOL_VERSION);
    serializedSoK << sigOfKnowledge;

    std::cout << "- Unserializing the Signature of Knowledge..." << endl;

    // unserialize the CDataStream object into a fresh SoK object
    SerialNumberSoK_small newSigOfKnowledge(ZCParams);
    serializedSoK >> newSigOfKnowledge;

    std::cout << "- Verifying the Signature of Knowledge..." << endl;

    // verify the signature of the received SoK
    bool res1 = newSigOfKnowledge.Verify(newCoin.getSerialNumber(), commitment.getCommitmentValue(), msghash);

    BOOST_CHECK_MESSAGE(res1,"SerialNumbwerSoK_small:: Verification failed\n");

    // if we made this far, all is good. Cheer the developer
    if(res1)
        std::cout << "  [YES] : [SoK is Ok! You did it champ. Have a cigar.]" << endl;

    std::cout << endl;
}

/*
BOOST_AUTO_TEST_CASE(signature_of_knowledge_tests)
{
    std::cout << "*** signature_of_knowledge_tests ***" << endl;
    std::cout << "------------------------------------" << endl;

    SelectParams(CBaseChainParams::MAIN);
    ZerocoinParams *ZCParams = Params().Zerocoin_Params(false);
    (void)ZCParams;

    // create message hash
    CHashWriter hasher(0,0);
    hasher << std::string("Chancellor on brink of second bailout for banks");
    uint256 msghash = hasher.GetHash();

    // mint a coin
    PrivateCoin newCoin(ZCParams, CoinDenomination::ZQ_ONE);

    // commit to this coin
    const CBigNum newCoin_value = newCoin.getPublicCoin().getValue();
    Commitment commitment(&(ZCParams->serialNumberSoKCommitmentGroup), newCoin_value);

    std::cout << "- Creating the Signature of Knowledge..." << endl;

    // create the signature of knowledge
    SerialNumberSoK_small sigOfKnowledge(ZCParams, newCoin, commitment, msghash);

    std::cout << "- Serializing the Signature of Knowledge..." << endl;

    // serialize the SoK to a CDataStream object.
    CDataStream serializedSoK(SER_NETWORK, PROTOCOL_VERSION);
    serializedSoK << sigOfKnowledge;

    std::cout << "- Unserializing the Signature of Knowledge..." << endl;

    // unserialize the CDataStream object into a fresh SoK object
    SerialNumberSoK_small newSigOfKnowledge(ZCParams);
    serializedSoK >> newSigOfKnowledge;

    std::cout << "- Verifying the Signature of Knowledge..." << endl;

    // verify the signature of the received SoK
    //bool res1 = newSigOfKnowledge.Verify(newCoin.getSerialNumber(), commitment.getCommitmentValue(), msghash);
    bool res1 = newSigOfKnowledge.Verify(newCoin.getSerialNumber(), commitment.getCommitmentValue(), msghash);

    BOOST_CHECK_MESSAGE(res1,"SerialNumbwerSoK_small:: Verification failed\n");

    std::cout << "- Checking that Verify fails for a different msghash..." << endl;

    // verify the signature on a random msghash
    CBigNum random_num = CBigNum::randBignum(CBigNum(2).pow(256));
    uint256 random_msghash = random_num.getuint256();

    bool res2 = newSigOfKnowledge.Verify(newCoin.getSerialNumber(), commitment.getCommitmentValue(), random_msghash);

    BOOST_CHECK_MESSAGE(!res2,"SerialNumbwerSoK_small:: Verification passed for random msghash\n");

    // if we made this far, all is good. Cheer the developer
    if(res1 && !res2)
        std::cout << "  [YES] : [SoK is Ok! You did it champ. Have a cigar.]" << endl;

    std::cout << endl;
}


BOOST_AUTO_TEST_CASE(signature_of_knowledge_benchmark)
{
    std::cout << "*** signature_of_knowledge_benchmark ***" << endl;
    std::cout << "----------------------------------------" << endl;

    SelectParams(CBaseChainParams::MAIN);
    ZerocoinParams *ZCParams = Params().Zerocoin_Params(false);
    (void)ZCParams;

    // create message hash
    CHashWriter hasher(0,0);
    hasher << std::string("Chancellor on brink of second bailout for banks");
    uint256 msghash = hasher.GetHash();

    // mint a coin
    PrivateCoin newCoin(ZCParams, CoinDenomination::ZQ_ONE);

    // commit to this coin
    const CBigNum newCoin_value = newCoin.getPublicCoin().getValue();
    Commitment commitment(&(ZCParams->serialNumberSoKCommitmentGroup), newCoin_value);

    // create the signature of knowledge (old)
    clock_t start_time = clock();
    SerialNumberSignatureOfKnowledge oldSigOfKnowledge(ZCParams, newCoin, commitment, msghash);
    clock_t total_time = clock() - start_time;
    double oldSok_time1 = total_time*1000.0/CLOCKS_PER_SEC;

    // create the signature of knowledge (new)
    start_time = clock();
    SerialNumberSoK_small sigOfKnowledge(ZCParams, newCoin, commitment, msghash);
    total_time = clock() - start_time;
    double newSok_time1 = total_time*1000.0/CLOCKS_PER_SEC;

    // verify the signature of knowledge (old)
    start_time = clock();
    oldSigOfKnowledge.Verify(newCoin.getSerialNumber(), commitment.getCommitmentValue(), msghash);
    total_time = clock() - start_time;
    double oldSok_time2 = total_time*1000.0/CLOCKS_PER_SEC;

    // verify the signature of knowledge (new)
    start_time = clock();
    sigOfKnowledge.Verify(newCoin.getSerialNumber(), commitment.getCommitmentValue(), msghash);
    total_time = clock() - start_time;
    double newSok_time2 = total_time*1000.0/CLOCKS_PER_SEC;

    // serialize the SoKs to CDataStream objects.
    CDataStream serializedSoK(SER_NETWORK, PROTOCOL_VERSION);
    serializedSoK << sigOfKnowledge;
    CDataStream serializedOldSoK(SER_NETWORK, PROTOCOL_VERSION);
    serializedOldSoK << oldSigOfKnowledge;
    std::cout << "- Size of the Serialized SoK (old): " << serializedOldSoK.size() << " bytes" << endl;
    std::cout << "- Time to sign (old): " << oldSok_time1 << " msec" << endl;
    std::cout << "- Time to verify (old): " << oldSok_time2 << " msec" << endl;
    std::cout << "- Size of the Serialized SoK (new): " << serializedSoK.size() << " bytes" << endl;
    std::cout << "- Time to sign (new): " << newSok_time1 << " msec" << endl;
    std::cout << "- Time to verify (new): " << newSok_time2 << " msec" << endl;

    std::cout << endl;
}


BOOST_AUTO_TEST_CASE(batch_signature_of_knowledge_tests)
{
    std::cout << "*** batch_signature_of_knowledge_tests ***" << endl;
    std::cout << "------------------------------------------" << endl;

    SelectParams(CBaseChainParams::MAIN);
    ZerocoinParams *ZCParams = Params().Zerocoin_Params(false);
    (void)ZCParams;

    // create message hash
    CHashWriter hasher(0,0);
    hasher << std::string("Chancellor on brink of second bailout for banks");
    uint256 msghash = hasher.GetHash();

    for(unsigned int k=1; k<=20; k++) {

        // mint k coins
        std::vector<PrivateCoin> coinlist;
        for(unsigned int i=0; i<k; i++) {
            PrivateCoin newCoin(ZCParams, CoinDenomination::ZQ_ONE);
            coinlist.push_back(newCoin);
        }

        // commit to these coins
        std::vector<Commitment> commitmentlist;
        for(unsigned int i=0; i<k; i++) {
            const CBigNum newCoin_value = coinlist[i].getPublicCoin().getValue();
            Commitment commitment(&(ZCParams->serialNumberSoKCommitmentGroup), newCoin_value);
            commitmentlist.push_back(commitment);
        }

        std::cout << "- Creating array of " << k << " Signatures of Knowledge..." << endl;

        // create k signatures of knowledge
        std::vector<SerialNumberSoK_small> siglist;

//clock_t start_time = clock();

        for(unsigned int i=0; i<k; i++) {
            SerialNumberSoK_small sigOfKnowledge(ZCParams, coinlist[i], commitmentlist[i], msghash);
            siglist.push_back(sigOfKnowledge);
        }

        std::cout << "- Packing and serializing the Signatures..." << endl;

        // pack the signature of knowledge
        std::vector<SerialNumberSoKProof> proofs;
        for(unsigned int i=0; i<k; i++) {
            SerialNumberSoKProof proof(siglist[i], coinlist[i].getSerialNumber(), commitmentlist[i].getCommitmentValue(), msghash);
            proofs.push_back(proof);
        }

        // serialize the proofs to a CDataStream object.
        std::vector<CDataStream> serializedProofs(proofs.size(), CDataStream(SER_NETWORK, PROTOCOL_VERSION));
        for(unsigned int i=0; i<serializedProofs.size(); i++) {
            serializedProofs[i] << proofs[i];
        }

        std::cout << "- Unserializing the Signature of Knowledge..." << endl;

        // unserialize the CDataStream object into a fresh SoK object
        std::vector<SerialNumberSoKProof> newproofs(serializedProofs.size(), SerialNumberSoKProof(ZCParams));
        for(unsigned int i=0; i<serializedProofs.size(); i++) {
            serializedProofs[i] >> newproofs[i];
        }


clock_t total_time = clock() - start_time;
double passed = total_time*1000.0/CLOCKS_PER_SEC;
std::cout << "----->" << passed << " msec" << endl;


        std::cout << "- Verifying the Signature of Knowledge..." << endl;

//start_time = clock();

        // verify the signature of the received SoK
        bool res1 = SerialNumberSoKProof::BatchVerify(newproofs);

        BOOST_CHECK_MESSAGE(res1,"SerialNumbwerSoK_small:: Verification failed\n");

        // if we made this far, all is good. Cheer the developer
        if(res1)
            std::cout << "  [YES] : [SoK is Ok! You did it champ. Have a cigar.]" << endl;

total_time = clock() - start_time;
passed = total_time*1000.0/CLOCKS_PER_SEC;
std::cout << "----->" << passed << " msec" << endl;

        std::cout << endl;

    }
}
*/
BOOST_AUTO_TEST_SUITE_END()
