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
        
        if( aParam.has( "manual-positioning" ) )
        {
            fManualPlacement = aParam["manual-positioning"]().as_bool();
        }
        if( aParam.has( "position-x" ) )
        {
            fCenterPositionX = aParam["position-x"]().as_double();
        }
        if( aParam.has( "position-y" ) )
        {
            fCenterPositionY = aParam["position-y"]().as_double();
        }
        if( aParam.has( "position-z" ) )
        {
            fCenterPositionZ = aParam["position-z"]().as_double();
        }
        if( aParam.has( "polarization-x" ) )
        {
            fPolarizationDirectionX = aParam["polarization-x"]().as_double();
        }
        if( aParam.has( "polarization-y" ) )
        {
            fPolarizationDirectionY = aParam["polarization-y"]().as_double();
        } 
        if( aParam.has( "polarization-z" ) )
        {
            fPolarizationDirectionZ = aParam["polarization-z"]().as_double();
        }
        if( aParam.has( "normal-x" ) )
        {
            fNormalDirectionX = aParam["normal-x"]().as_double();
        } 
        if( aParam.has( "normal-y" ) )
        {
            fNormalDirectionY = aParam["normal-y"]().as_double();
        } 
        if( aParam.has( "normal-z" ) )
        {
            fNormalDirectionZ = aParam["normal-z"]().as_double();
        } 
        
        if(fPolarizationDirectionX==fPolarizationDirectionY
			&& fPolarizationDirectionX==fPolarizationDirectionZ
			&& fPolarizationDirectionX==0.0)
		{
            fPolarizationDirectionY=-1.0;
        } 
           
        if(fNormalDirectionX==fNormalDirectionY==fNormalDirectionZ==0.0)
        {
            if(fCenterPositionX==fCenterPositionY==fCenterPositionZ==0.0) 
            {
                fNormalDirectionX = 1.0;
            } else
            {
                double distance = sqrt(fCenterPositionX*fCenterPositionX + 
                                        fCenterPositionY*fCenterPositionY +
                                        fCenterPositionZ*fCenterPositionZ);
                                        
                fNormalDirectionX = -fCenterPositionX/distance;
                fNormalDirectionY = -fCenterPositionY/distance;
                fNormalDirectionZ = -fCenterPositionZ/distance;
            }
        }
        
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
        if(fManualPlacement)
            PlaceAtConfigured(modelElement);
        else
            PlaceAtPosition(modelElement, elementRadius, theta, zPosition);
	}
    
    void AntennaElementPositioner::PlaceAtPosition(Receiver &modelElement, double elementRadius, double theta, double zPosition)
	{
		modelElement.SetCenterPosition({elementRadius * cos(theta) , elementRadius * sin(theta) , zPosition });
		modelElement.SetPolarizationDirection({sin(theta), -cos(theta), 0.0});
		modelElement.SetCrossPolarizationDirection({0.0, 0.0, 1.0});  // longitudinal axis of array.
		modelElement.SetNormalDirection({-cos(theta), -sin(theta), 0.0}); //Say normals point inwards
	}

	void AntennaElementPositioner::PlaceAtConfigured(Receiver &modelElement)
	{
        LDEBUG(lmclog, "Placing element manually");
        
        modelElement.SetCenterPosition({fCenterPositionX , fCenterPositionY , fCenterPositionZ });
		modelElement.SetPolarizationDirection({fPolarizationDirectionX, fPolarizationDirectionY, fPolarizationDirectionZ});
		modelElement.SetCrossPolarizationDirection({0.0, 0.0, 1.0});  // longitudinal axis of array.
		modelElement.SetNormalDirection({fNormalDirectionX, fNormalDirectionY, fNormalDirectionZ}); //Say normals point inwards
	}

} /* namespace locust */

