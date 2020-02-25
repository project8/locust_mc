/*
 * LMCDipoleAntenna.cc
 *
 *  Created on: Jan 24, 2020
 *      Author: pslocum
 */

#include "LMCDipoleAntenna.hh"
#include "logger.hh"
using std::string;



namespace locust
{
    LOGGER( lmclog, "DipoleAntenna" );

    DipoleAntenna::DipoleAntenna():
    fMomentVector( 0., 0., 1.0 ),
    fMagneticDipole( true )
    {
    }

    DipoleAntenna::~DipoleAntenna()
    {
    }

    bool DipoleAntenna::Configure( const scarab::param_node& aParam )
    {

	if( !TransmitterHardware::Configure(aParam))
	{
     	    LERROR(lmclog,"Error configuring TransmitterHardware class from DipoleAntenna child class");
	}

        if( aParam.has( "dipoleantenna-momentX" ) )
        {
            fMomentVector.SetX(aParam["dipoleantenna-momentX"]().as_double());
        }

        if( aParam.has( "dipoleantenna-momentY" ) )
        {
            fMomentVector.SetY(aParam["dipoleantenna-momentY"]().as_double());
        }

        if( aParam.has( "dipoleantenna-momentZ" ) )
        {
            fMomentVector.SetZ(aParam["dipoleantenna-momentZ"]().as_double());
        }

        if( aParam.has( "dipoleantenna-magnetic" ) )
        {
            fMagneticDipole = aParam["dipoleantenna-magnetic"]().as_bool();
        }

    	return true;
    }

    void DipoleAntenna::TxHardwareSayHello()
     {
    	printf("Dipole says hello\n"); getchar();
     }

    double DipoleAntenna::GetPatternFactor(LMCThreeVector pointOfInterest, int antennaNumber)
    {
    	double patternFactor=0.0; 
    	if (fMagneticDipole)
    	{
		patternFactor = pointOfInterest.Unit().Cross(fMomentVector.Unit()).Magnitude(); // sin(theta)
    	}
    	else
    	{
		patternFactor = pointOfInterest.Unit().Cross(fMomentVector.Unit()).Magnitude(); // sin(theta)
    	}
    	return patternFactor;

    }





} /* namespace locust */
