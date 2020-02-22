/*
 * LMCTurnstileAntenna.cc
 *
 *  Created on: Jan 24, 2020
 *      Author: pslocum
 */

#include "LMCTurnstileAntenna.hh"
#include "logger.hh"
using std::string;



namespace locust
{

	LOGGER( lmclog, "TurnstileAntenna" );

    TurnstileAntenna::TurnstileAntenna()
    {
    }

    TurnstileAntenna::~TurnstileAntenna()
    {
    }

    bool TurnstileAntenna::Initialize()
    {

    	SetNAntennas( 2 );
    	fMomentVector.resize(GetNAntennas());

    	fMomentVector[0].SetX(1.);
    	fMomentVector[0].SetY(0.);
    	fMomentVector[0].SetZ(0.);

    	fMomentVector[1].SetX(0.);
    	fMomentVector[1].SetY(1.);
    	fMomentVector[1].SetZ(0.);

    	return true;
    }

    bool TurnstileAntenna::Configure( const scarab::param_node& aParam )
    {
    	Initialize();
    	return true;
    }

    void TurnstileAntenna::TxHardwareSayHello()
     {
    	printf("Turnstile says hello\n"); getchar();
     }

    double TurnstileAntenna::GetPatternFactor(LMCThreeVector pointOfInterest, int antennaNumber)
    {
    	// sin(theta)
    	double patternFactor = pointOfInterest.Unit().Cross(fMomentVector[antennaNumber].Unit()).Magnitude();
    	return patternFactor;

    }





} /* namespace locust */
