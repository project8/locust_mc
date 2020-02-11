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



    double KassTransmitter::GetAOIFactor(LMCThreeVector IncidentKVector, double ElementPhi)
    {
        LMCThreeVector ElementNormalVector;
        ElementNormalVector.SetComponents(cos(ElementPhi), sin(ElementPhi), 0.0);
        double AOIFactor = fabs(IncidentKVector.Unit().Dot(ElementNormalVector));
        //printf("cos aoi is %f\n", AOIFactor);
        return AOIFactor;
    }


    // EField cross pol with aoi dot product, at element.
    double KassTransmitter::GetEFieldCoPol(Receiver* currentElement, LMCThreeVector IncidentElectricField, LMCThreeVector IncidentKVector, double ElementPhi)
    {
        double AOIFactor = GetAOIFactor(IncidentKVector, ElementPhi);  // k dot elementnormal
        LMCThreeVector ElementPolarizationVector = currentElement->GetPolarizationDirection();
        double EFieldCoPol = IncidentElectricField.Dot(ElementPolarizationVector) * AOIFactor;

        return EFieldCoPol;
    }


    double KassTransmitter::GetEFieldCrossPol(Receiver* currentElement, LMCThreeVector IncidentElectricField, LMCThreeVector IncidentKVector, double ElementPhi)
    {
        double AOIFactor = GetAOIFactor(IncidentKVector, ElementPhi);  // k dot elementnormal
        LMCThreeVector ElementCrossPolarizationVector;
        ElementCrossPolarizationVector.SetComponents(0.0, 0.0, 1.0);  // axial direction.
        double EFieldCrossPol = IncidentElectricField.Dot(ElementCrossPolarizationVector) * AOIFactor;

        return EFieldCrossPol;
    }

    void KassTransmitter::InitializeFieldPoints(std::vector< Channel<Receiver*> > allRxChannels)
    {

    	int nChannels = allRxChannels.size();
    	int nReceivers = allRxChannels[0].size();

        for(int channelIndex = 0; channelIndex < nChannels; ++channelIndex)
        {
            for(int elementIndex = 0; elementIndex < nReceivers; ++elementIndex)
            {
            	Receiver* currentElement = allRxChannels[channelIndex][elementIndex];
            	fFieldSolver.AddFieldPoint(currentElement->GetPosition());
            }
        }

    }


	double* KassTransmitter::SolveKassFields(Receiver* currentElement, double ElementPhi, double tReceiverTime, unsigned tTotalElementIndex)
    {

        fFieldSolver.SetFieldEvent(tReceiverTime, tTotalElementIndex);
        fFieldSolver.SolveFieldSolutions();

        LMCThreeVector tRadiatedElectricField = fFieldSolver.GetElectricField();
        LMCThreeVector tRadiatedMagneticField = fFieldSolver.GetMagneticField();
        locust::Particle tCurrentParticle = fFieldSolver.GetRetardedParticle();

        LMCThreeVector tDirection = currentElement->GetPosition() - tCurrentParticle.GetPosition(true);
	    double tVelZ = tCurrentParticle.GetVelocity(true).Z();
	    double tCosTheta =  tVelZ * tDirection.Z() /  tDirection.Magnitude() / fabs(tVelZ);
	    double tDopplerFrequency  = tCurrentParticle.GetCyclotronFrequency() / ( 1. - fabs(tVelZ) / LMCConst::C() * tCosTheta);

	    double tEFieldCoPol = GetEFieldCoPol(currentElement, tRadiatedElectricField, tRadiatedElectricField.Cross(tRadiatedMagneticField), ElementPhi);
	    double tEFieldCrossPol = GetEFieldCrossPol(currentElement, tRadiatedElectricField, tRadiatedElectricField.Cross(tRadiatedMagneticField), ElementPhi);

        double* tSolution = new double[2];
        tSolution[0] = tEFieldCoPol;
        tSolution[1] = tDopplerFrequency;

	    return tSolution;

    }



} /* namespace locust */
