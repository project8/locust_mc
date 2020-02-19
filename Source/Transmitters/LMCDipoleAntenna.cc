/*
 * LMCDipoleAntenna.cc
 *
 *  Created on: Jan 24, 2020
 *      Author: pslocum
 */

#include "LMCDipoleAntenna.hh"


namespace locust
{

    DipoleAntenna::DipoleAntenna():
    fMomentVector( 0., 0., 1.0 )
    {
    }

    DipoleAntenna::~DipoleAntenna()
    {
    }


    void DipoleAntenna::TxHardwareSayHello()
     {
    	printf("Dipole says hello\n"); getchar();
     }

    double DipoleAntenna::GetPatternFactor(LMCThreeVector pointOfInterest)
    {
    	return pointOfInterest.Unit().Dot(fMomentVector.Orthogonal().Unit()); // cos(theta) dependence
    }





} /* namespace locust */
