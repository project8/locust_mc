/*
 * LMCSlottedWaveguide.cc
 *
 *  Created on: Feb 25, 2020
 *      Author: pslocum
 */

#include "LMCSlottedWaveguide.hh"
using std::string;


namespace locust
{

	LOGGER( lmclog, "SlottedWaveguide" );

    SlottedWaveguide::SlottedWaveguide():
    		fImpedanceTransformation( 0.339 ) // sqrt(50./435.) unless already included in HFSS TF.
    {
    }

    SlottedWaveguide::~SlottedWaveguide()
    {
    }

    bool SlottedWaveguide::Configure( const scarab::param_node& aParam )
    {

    	if( !PowerCombiner::Configure(aParam))
    	{
    		LERROR(lmclog,"Error configuring PowerCombiner class from SlottedWaveguide subclass");
    		return false;
    	}

        if( aParam.has( "impedance-transformation" ) )
        {
            fImpedanceTransformation = aParam["impedance-transformation"]().as_double();
        }

    	SetVoltageDampingFactors();
    	return true;
    }



	bool SlottedWaveguide::SetVoltageDampingFactors()
	{

		for (unsigned z_index=0; z_index<GetNElementsPerStrip(); z_index++)
		{
			double aFactor = fImpedanceTransformation;
			//Adhoc scaling in case the slotted waveguide has fewer/more slots than 10
			aFactor = aFactor/sqrt(GetNElementsPerStrip()/10.);
			SetDampingFactor(z_index, aFactor);
		}

		return true;

	}

    Receiver* SlottedWaveguide::ChooseElement()
    {
    	SlotAntenna* aSlot = new SlotAntenna;
    	return aSlot;
    }



} /* namespace locust */
