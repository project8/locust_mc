/*
 * LMCPlaneWaveTransmitter.cc
 *
 *  Created on: May 4, 2019
 *      Author: pslocum
 */

#include "LMCPlaneWaveTransmitter.hh"
#include "logger.hh"

using std::string;


namespace locust
{
    LOGGER( lmclog, "PlaneWaveTransmitter" );

    PlaneWaveTransmitter::PlaneWaveTransmitter():
		fAOI( 0.),
		fAmplitude( 0.0 ),
		fRF_Frequency( 0. )
    {
    }

    PlaneWaveTransmitter::~PlaneWaveTransmitter()
    {
    }

    bool PlaneWaveTransmitter::Configure( const scarab::param_node& aParam )
    {

		if( aParam.has( "transmitter-frequency" ) )
		{
            fRF_Frequency= aParam["transmitter-frequency"]().as_double();
		}
		if( aParam.has( "AOI" ) )
		{
            fAOI= aParam["AOI"]().as_double()*2.*LMCConst::Pi()/360.;

		}
		if( aParam.has( "planewave-amplitude" ) )
		{
            fAmplitude= aParam["planewave-amplitude"]().as_double();
		}


        return true;
    }

    void PlaneWaveTransmitter::AddPropagationPhaseDelay(LMCThreeVector pointOfInterest)
    {
	//Assuming the element strip is always along Z
      
	/*double meanZ = GetMeanofFieldPoints(2);
	  double distanceFromCenter = meanZ - pointOfInterest.GetZ();*/
      
	/*int z_index = fieldPointIndex%nElementsPerStrip;
    	double stripLength = (nElementsPerStrip-1)*elementSpacing;
    	double distanceFromCenter = stripLength/2. - z_index*elementSpacing;*/

      double endOfStrip = GetFieldPoint(0).GetZ();
      double distanceFromCenter = pointOfInterest.GetZ() - endOfStrip;
    
    	double phaseDelay = 2*LMCConst::Pi()*distanceFromCenter*sin(fAOI)*fRF_Frequency/LMCConst::C();
	Transmitter::AddPropagationPhaseDelay(phaseDelay);
    }

    void PlaneWaveTransmitter::AddIncidentKVector(LMCThreeVector pointOfInterest)
    {
	LMCThreeVector incidentKVector(cos(fAOI), 0.0, sin(fAOI));
	Transmitter::AddIncidentKVector(incidentKVector);
    	//fIncidentKVector.SetComponents(cos(fAOI), 0.0, sin(fAOI));
    }

    double* PlaneWaveTransmitter::GetEFieldCoPol(int fieldPointIndex, double dt)
    {
    	double initialPhaseDelay = GetPropagationPhaseDelay(fieldPointIndex); 
		double fieldAmp = fAmplitude;

		if (fieldPointIndex==0) fPhaseDelay += 2. * LMCConst::Pi() * fRF_Frequency * dt;
		//if ( (zIndex == fieldPointIndex0) && (channelIndex == 0) ) fPhaseDelay += 2. * LMCConst::Pi() * fRF_Frequency * dt;
		double fieldValue = fieldAmp*cos(fPhaseDelay + initialPhaseDelay);
		//AddIncidentKVector(pointOfInterest);

        double* fieldSolution = new double[2];
        fieldSolution[0] = fieldValue;
        fieldSolution[1] = 2. * LMCConst::Pi() * fRF_Frequency;  // rad/s

        return fieldSolution;
    }



} /* namespace locust */
