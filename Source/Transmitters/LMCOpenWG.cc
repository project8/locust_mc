/*
 * LMCOpenWG.cc
 *
 *  Created on: Apr 6, 2022
 *      Author: P. T. Surukuchi
 * Needs implementation still 
 */

#include "LMCOpenWG.hh"
#include "logger.hh"
using std::string;



namespace locust
{
    LOGGER( lmclog, "OpenWG" );

    OpenWG::OpenWG():
    fMomentVector( 0., 0., 1.0 )
    {
    }

    OpenWG::~OpenWG()
    {
    }

    bool OpenWG::Configure( const scarab::param_node& aParam )
    {

	if( !TransmitterHardware::Configure(aParam))
	{
     	    LERROR(lmclog,"Error configuring TransmitterHardware class from OpenWG child class");
	}

        if( aParam.has( "OpenWG-momentX" ) )
        {
            fMomentVector.SetX(aParam["OpenWG-momentX"]().as_double());
        }

        if( aParam.has( "OpenWG-momentY" ) )
        {
            fMomentVector.SetY(aParam["OpenWG-momentY"]().as_double());
        }

        if( aParam.has( "OpenWG-momentZ" ) )
        {
            fMomentVector.SetZ(aParam["OpenWG-momentZ"]().as_double());
        }

    	return true;
    }

    void OpenWG::TxHardwareSayHello()
     {
    	printf("Dipole says hello\n"); getchar();
     }

    double OpenWG::GetPatternFactor(LMCThreeVector pointOfInterest, int antennaNumber)
    {
    	double patternFactor=0.0; 
	LMCThreeVector incidentKVector=ExtractIncidentKVector(pointOfInterest);
    //sin(theta) between the line joining the dipole to the POI and the moment
    double patternFactor_phi = incidentKVector.Unit().Cross(fMomentVector.Unit()).Magnitude(); // sin(theta)
    //The Efield is in phi direction w.r.t the dipole. 
    //To calculate the copol, the component of the phi in copol direction has to be taken
    //This is done by  projecting the moment of the dipole in z direction of the array
    //Fails miserably if the copol of the antenna/slots is not in the standard direction
    patternFactor=patternFactor_phi*LMCThreeVector(0,0,1).Dot(fMomentVector.Unit());
    return patternFactor;

    }





} /* namespace locust */
