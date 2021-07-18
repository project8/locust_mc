/*
 * LMCCavityModes.cc
 *
 *  Created on: Jul 16, 2021
 *      Author: pslocum
 */

#include "LMCCavityModes.hh"
using std::string;


namespace locust
{

	LOGGER( lmclog, "CavityModes" );

    CavityModes::CavityModes()
    {
    }

    CavityModes::~CavityModes()
    {
    }

    bool CavityModes::SetProbeLocations()
    {
    	std::vector<double> probeZ;
    	probeZ.resize(1);
    	probeZ[0] = 0.;  // TO-DO:  parametrize this.
    	SetCavityProbeZ(probeZ);

    	std::vector<double> probeTheta;
    	probeTheta.resize(1);
    	probeTheta[0] = 0.;  // TO-DO:  parametrize this.
    	SetCavityProbeTheta(probeTheta);

    	return true;
    }

    bool CavityModes::Configure( const scarab::param_node& aParam )
    {

    	if( !PowerCombiner::Configure(aParam))
    	{
    		LERROR(lmclog,"Error configuring PowerCombiner class from CavityModes subclass");
    		return false;
    	}

    	if ( aParam.has( "n-cavity-probes" ) )
    	{
    		SetNCavityProbes(aParam["n-cavity-probes"]().as_int());
    	}

    	if ( aParam.has( "cavity-probe-impedance" ) )
    	{
    		SetCavityProbeImpedance(aParam["cavity-probe-impedance"]().as_double());
    	}

    	SetProbeLocations();

    	return true;
    }



} /* namespace locust */
