/*
 * LMCAntennaElementPositioner.cc
 *
 *  Created on: Jul 22, 2020
 *      Author: pslocum
 */

#include "LMCAntennaElementPositioner.hh"
#include <iostream>
#include "logger.hh"

namespace locust
{

	LOGGER( lmclog, "AntennaElementPositioner" );


    AntennaElementPositioner::AntennaElementPositioner() {}
    AntennaElementPositioner::~AntennaElementPositioner() {}


    bool AntennaElementPositioner::Configure( const scarab::param_node& aParam )
    {
        return true;
    }

    double AntennaElementPositioner::GetPositionZ(double zShiftArray, int channelIndex, int nChannels,
    		int nSubarrays, int nReceivers, double elementSpacingZ, int receiverIndex)
    {
		double zPosition =  zShiftArray +
				(int(channelIndex/(nChannels/nSubarrays))-((nSubarrays -1.)/2.) )*nReceivers*elementSpacingZ +
				(receiverIndex - (nReceivers - 1.) /2.) * elementSpacingZ;
    	return zPosition;
    }

    double AntennaElementPositioner::GetTheta(int channelIndex, double dThetaArray)
    {
    	double tTheta = channelIndex * dThetaArray;
    	return tTheta;
    }

	void AntennaElementPositioner::PlaceElement(Receiver &modelElement, double elementRadius, double theta, double zPosition)
	{
		modelElement.SetCenterPosition({elementRadius * cos(theta) , elementRadius * sin(theta) , zPosition });
		modelElement.SetPolarizationDirection({sin(theta), -cos(theta), 0.0});
		modelElement.SetCrossPolarizationDirection({0.0, 0.0, 1.0});  // longitudinal axis of array.
		modelElement.SetNormalDirection({-cos(theta), -sin(theta), 0.0}); //Say normals point inwards
	}



} /* namespace locust */

