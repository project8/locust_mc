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
	LMCThreeVector incidentKVector=ExtractIncidentKVector(pointOfInterest);
    	if (fMagneticDipole)
    	{
		//sin(theta) between the line joining the dipole to the POI and the moment
		double patternFactor_phi = incidentKVector.Unit().Cross(fMomentVector.Unit()).Magnitude(); // sin(theta)
		//The Efield is in phi direction w.r.t the dipole. 
		//To calculate the copol, the component of the phi in copol direction has to be taken
		//This is done by  projecting the moment of the dipole in z direction of the array
		//Fails miserably if the copol of the antenna/slots is not in the standard direction
		patternFactor=patternFactor_phi*LMCThreeVector(0,0,1).Dot(fMomentVector.Unit());
    	}
    	else //If electric dipole
    	{
		//sin(theta) between the line joining the dipole to the POI and the moment
		//The actual factor is cos((π/2)cos(θ))/sin(θ), should be implemented later
		double patternFactor_theta = incidentKVector.Unit().Cross(fMomentVector.Unit()).Magnitude(); // sin(theta)
		//The Efield is in theta direction w.r.t the dipole. 
		//To calculate the copol, the component of the phi in copol direction has to be taken
		//This is done by  taking the sin of angle between the moment and the z direction
		//Fails miserably if the copol of the antenna/slots is not in the standard direction
		patternFactor=patternFactor_theta*LMCThreeVector(0,0,1).Cross(fMomentVector.Unit()).Magnitude();
    	}
    	return patternFactor;

    }





} /* namespace locust */
