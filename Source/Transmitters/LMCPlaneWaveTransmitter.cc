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


    double PlaneWaveTransmitter::GetAOIFactor(Receiver* currentElement)
    {
    	LMCThreeVector incidentKVector;
    	incidentKVector.SetComponents(cos(fAOI), 0.0, sin(fAOI));
    	return currentElement->GetPatternFactor(incidentKVector, *currentElement);
    }


    double PlaneWaveTransmitter::GetPWPhaseDelayAtPatch(int z_index, double elementSpacing, int nElementsPerStrip)
    {

    	double stripLength = (nElementsPerStrip-1)*elementSpacing;
    	double distanceFromCenter = stripLength/2. - z_index*elementSpacing;

    	double phasedelay = 2*LMCConst::Pi()*distanceFromCenter*sin(fAOI)*fRF_Frequency/LMCConst::C();

    	return phasedelay;
    }



    double* PlaneWaveTransmitter::GetEFieldCoPol(Receiver* currentElement, int channelIndex, int zIndex, double elementSpacing, int nElementsPerStrip, double dt)
    {

    	double initialPhaseDelay = GetPWPhaseDelayAtPatch(zIndex, elementSpacing, nElementsPerStrip);
		double fieldAmp = fAmplitude*GetAOIFactor(currentElement);

		if ( (zIndex == 0) && (channelIndex == 0) ) fPhaseDelay += 2. * LMCConst::Pi() * fRF_Frequency * dt;
		double fieldValue = fieldAmp*cos(fPhaseDelay + initialPhaseDelay);

        double* fieldSolution = new double[2];
        fieldSolution[0] = fieldValue;
        fieldSolution[1] = 2. * LMCConst::Pi() * fRF_Frequency;  // rad/s

        return fieldSolution;
    }



} /* namespace locust */
