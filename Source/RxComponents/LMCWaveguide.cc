/*
 * LMCWaveguide.cc
 *
 *  Created on: Nov 15, 2021
 *      Author: atelles
 */

#include "LMCWaveguide.hh"
using std::string;


namespace locust
{

	LOGGER( lmclog, "Waveguide" );

    Waveguide::Waveguide():
    		fImpedanceTransformation( 0.263 ) // default sqrt(50./725.), where 725 ohms is the waveguide TF impedance
    {
    }

    Waveguide::~Waveguide()
    {
    }

    bool Waveguide::Configure( const scarab::param_node& aParam )
    {

    	if( !PowerCombiner::Configure(aParam))
    	{
    		LERROR(lmclog,"Error configuring PowerCombiner class from Waveguide subclass");
    		return false;
    	}

        if( aParam.has( "impedance-transformation" ) )
        {
            fImpedanceTransformation = aParam["impedance-transformation"]().as_double();
        }

    	return true;
    }

} /* namespace locust */
