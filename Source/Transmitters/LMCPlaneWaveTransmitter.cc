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
    	LMCThreeVector IncidentKVector;
    	IncidentKVector.SetComponents(cos(fAOI), 0.0, sin(fAOI));
    	double AOIFactor = fabs(IncidentKVector.Dot(currentElement->GetNormalDirection()));
    	return AOIFactor;
    }


    double PlaneWaveTransmitter::GetPWPhaseDelayAtPatch(int z_index, double elementSpacing, int nElementsPerStrip)
    {
    	double phasedelay = 0.;
    	if(fAOI >= 0)
    	{
    		phasedelay = 2*LMCConst::Pi()*z_index*elementSpacing*sin(fAOI)*fRF_Frequency/LMCConst::C();
    	}
    	else
    	{
    		phasedelay = (nElementsPerStrip - z_index)*2*LMCConst::Pi()*elementSpacing*sin(fAOI)*fRF_Frequency/LMCConst::C();
    	}
    	return phasedelay;
    }



    double* PlaneWaveTransmitter::GetEFieldCoPol(Receiver* currentElement, int z_index, double elementSpacing, int nElementsPerStrip, double dt)
    {

    	double initialPhaseDelay = GetPWPhaseDelayAtPatch(z_index, elementSpacing, nElementsPerStrip);
		double fieldAmp = fAmplitude*GetAOIFactor(currentElement);

		initialPhaseDelay = 0.;
		fieldAmp = fAmplitude;

		if (z_index == 0) fPhaseDelay += 2. * LMCConst::Pi() * fRF_Frequency * dt;
		double fieldValue = fieldAmp*cos(fPhaseDelay + initialPhaseDelay);


        double* fieldSolution = new double[2];
        fieldSolution[0] = fieldValue;
        fieldSolution[1] = 2. * LMCConst::Pi() * fRF_Frequency;  // rad/s

        return fieldSolution;
    }



} /* namespace locust */
