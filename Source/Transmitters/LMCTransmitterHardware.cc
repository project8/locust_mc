/*
 * LMCTransmitterHardware.cc
 *
 *  Created on: Feb 18, 2020
 *      Author: pslocum
 */

#include "LMCTransmitterHardware.hh"


namespace locust
{
    TransmitterHardware::TransmitterHardware():
    	fDrivePhaseDifference( 0. ),
    	fNAntennas( 1 ),
    	fAntennaPositionX( 0.0 ),
    	fAntennaPositionY( 0.0 ),
    	fAntennaPositionZ( 0.0 )
    {
    }

    TransmitterHardware::~TransmitterHardware() {}

    bool TransmitterHardware::Configure( const scarab::param_node& aParam )
    {

        if( aParam.has( "antenna-x-position" ) )
        {
            fAntennaPositionX= aParam["antenna-x-position"]().as_double();
        }
        
        if( aParam.has( "antenna-y-position" ) )
        {
            fAntennaPositionY = aParam["antenna-y-position"]().as_double();
        }
        
        if( aParam.has( "antenna-z-position" ) )
        {
            fAntennaPositionZ = aParam["antenna-z-position"]().as_double();
        }
        fAntennaPosition.SetComponents(fAntennaPositionX,fAntennaPositionY,fAntennaPositionZ);
        return true;
    }
        
    void TransmitterHardware::TxHardwareSayHello()
    {
    	printf("TransmitterHardware says hello\n"); getchar();
    }

    int TransmitterHardware::GetNAntennas()
    {
    	return fNAntennas;
    }

    void TransmitterHardware::SetNAntennas(int aNumber)
    {
    	fNAntennas = aNumber;
    }

    void TransmitterHardware::SetAntennaPosition(const LMCThreeVector &antennaPosition)
    {
        fAntennaPosition=antennaPosition;
    }

    LMCThreeVector TransmitterHardware::GetAntennaPosition() const
    {
        return fAntennaPosition;
    }
    
    LMCThreeVector TransmitterHardware::ExtractIncidentKVector(LMCThreeVector pointOfInterest)
    {
    	LMCThreeVector incidentKVector;

    	double relativeElementPosX=pointOfInterest.GetX() - fAntennaPosition.GetX();
        double relativeElementPosY=pointOfInterest.GetY() - fAntennaPosition.GetY();
        double relativeElementPosZ=pointOfInterest.GetZ() - fAntennaPosition.GetZ();
     	incidentKVector.SetComponents(relativeElementPosX, relativeElementPosY, relativeElementPosZ);
	return incidentKVector;
    }
    
    double TransmitterHardware::GetPropagationDistance(LMCThreeVector pointOfInterest)
    {
	return ExtractIncidentKVector(pointOfInterest).Magnitude();
    }

} /* namespace locust */

