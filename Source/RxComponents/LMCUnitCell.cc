/*
 * LMCUnitCell.cc
 *
 *  Created on: Feb 25, 2020
 *      Author: pslocum
 */

#include "LMCUnitCell.hh"
using std::string;


namespace locust
{

	LOGGER( lmclog, "UnitCell" );

    UnitCell::UnitCell()
    {
    }

    UnitCell::~UnitCell()
    {
    }

    bool UnitCell::Configure( const scarab::param_node& aParam )
    {

    	if( !PowerCombiner::Configure(aParam))
    	{
    		LERROR(lmclog,"Error configuring PowerCombiner class from UnitCell subclass");
    		return false;
    	}


    	if (GetNElementsPerStrip()%2 != 0)
    	{
    		LERROR(lmclog,"Center-fed power combining config expects even number of patches per strip.");
    		return false;
    	}


    	if(aParam["power-combining-feed"]().as_string() == "unit-cell-one-quarter")
    	{
			SetJunctionLoss(0.87);
			SetPatchLoss(0.38);
			SetAmplifierLoss(0.66);
			SetEndPatchLoss(0.95);
    	}

    	if(aParam["power-combining-feed"]().as_string() == "unit-cell-seven-eighths")
    	{
			SetJunctionLoss(0.75);
			SetPatchLoss(0.6);
			SetAmplifierLoss(0.66);
			SetEndPatchLoss(0.95);
    	}

    	if(aParam["power-combining-feed"]().as_string() == "unit-cell-nine-sixteenths")
    	{
			SetJunctionLoss(0.8);
			SetPatchLoss(0.52);
			SetAmplifierLoss(0.66);
			SetEndPatchLoss(0.95);
    	}

    	SetVoltageDampingFactors();
    	return true;
    }


	bool UnitCell::SetVoltageDampingFactors()
	{

		for (unsigned z_index=0; z_index<GetNElementsPerStrip(); z_index++)
		{
			int njunctions = fabs((double)z_index - GetNElementsPerStrip()/2.) - 1;
			if (z_index >= GetNElementsPerStrip()/2) njunctions += 1; // compensate for patches to the right of amp.
			if (z_index == 0 || z_index == GetNElementsPerStrip()-1)
	      	{
				double aFactor = GetEndPatchLoss()*pow(GetJunctionLoss(), njunctions)*GetAmplifierLoss();
				SetDampingFactor(z_index, aFactor);
	      	}
			else
	      	{
				double aFactor = GetPatchLoss()*pow(GetJunctionLoss(), njunctions)*GetAmplifierLoss();
				SetDampingFactor(z_index, aFactor);
	      	}
		}

		return true;

	}




} /* namespace locust */
