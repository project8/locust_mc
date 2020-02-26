/*
 * LMCKassTransmitter.cc
 *
 *  Created on: May 4, 2019
 *      Author: pslocum
 */

#include "LMCKassTransmitter.hh"
#include "logger.hh"
#include <thread>
#include <algorithm>

#include <iostream>
#include <fstream>
#include <math.h>

using std::string;


namespace locust
{
    LOGGER( lmclog, "KassTransmitter" );

    KassTransmitter::KassTransmitter()
    {
    }

    KassTransmitter::~KassTransmitter()
    {
    }

    bool KassTransmitter::Configure( const scarab::param_node& aParam )
    {
        return true;
    }

    bool KassTransmitter::IsKassiopeia()
    {
    	return true;
    }

    void KassTransmitter::InitializeFieldPoint(LMCThreeVector fieldPoint)
    {
	Transmitter::InitializeFieldPoint(fieldPoint);
    	fFieldSolver.AddFieldPoint(fieldPoint);
    }

    double* KassTransmitter::SolveKassFields(LMCThreeVector pointOfInterest, LMCThreeVector coPolDirection, double tReceiverTime, unsigned tTotalElementIndex)
    {

        fFieldSolver.SetFieldEvent(tReceiverTime, tTotalElementIndex);
        fFieldSolver.SolveFieldSolutions();

        LMCThreeVector tRadiatedElectricField = fFieldSolver.GetElectricField();
        LMCThreeVector tRadiatedMagneticField = fFieldSolver.GetMagneticField();
        locust::Particle tCurrentParticle = fFieldSolver.GetRetardedParticle();

        LMCThreeVector tDirection = pointOfInterest - tCurrentParticle.GetPosition(true);
	    double tVelZ = tCurrentParticle.GetVelocity(true).Z();
	    double tCosTheta =  tVelZ * tDirection.Z() /  tDirection.Magnitude() / fabs(tVelZ);
	    double tDopplerFrequency  = tCurrentParticle.GetCyclotronFrequency() / ( 1. - fabs(tVelZ) / LMCConst::C() * tCosTheta);

	    double tEFieldCoPol = tRadiatedElectricField.Dot(coPolDirection);
	    SetIncidentKVector(tTotalElementIndex,tRadiatedElectricField.Cross(tRadiatedMagneticField));

        double* tSolution = new double[2];
        tSolution[0] = tEFieldCoPol;
        tSolution[1] = tDopplerFrequency;

	    return tSolution;

    }



} /* namespace locust */
