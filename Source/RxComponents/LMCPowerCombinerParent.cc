/*
 * LMCPowerCombiner.cc
 *
 *  Created on: Feb 25, 2020
 *      Author: pslocum
 */

#include "LMCPowerCombinerParent.hh"
#include <iostream>
#include "logger.hh"

namespace locust
{

	LOGGER( lmclog, "PowerCombinerParent" );


    PowerCombinerParent::PowerCombinerParent():
			fnElementsPerStrip( 2 ),
			fdampingFactors( 0 ),
            fjunctionLoss( 1.0 ),
            fpatchLoss( 0.6 ),
            famplifierLoss( 0.66 ),
            fendPatchLoss( 1.0 ),
            fjunctionResistance( 0.3 )

    {}
    PowerCombinerParent::~PowerCombinerParent() {}


    bool PowerCombinerParent::Configure( const scarab::param_node& aParam )
    {
        if( aParam.has( "nelements-per-strip" ) )
        {
            fnElementsPerStrip = aParam["nelements-per-strip"]().as_int();
            fdampingFactors.resize( fnElementsPerStrip );
        }

        return true;

    }


    void PowerCombinerParent::SayHello()
    {
    	printf("powercombiner says hello\n"); getchar();
    }


	bool PowerCombinerParent::AddOneVoltageToStripSum(Signal* aSignal, double VoltageFIRSample, double phi_LO, unsigned z_index, unsigned sampleIndex)
	{

		VoltageFIRSample *= GetDampingFactor(z_index);
		aSignal->LongSignalTimeComplex()[sampleIndex][0] += 2.*VoltageFIRSample * sin(phi_LO);
		aSignal->LongSignalTimeComplex()[sampleIndex][1] += 2.*VoltageFIRSample * cos(phi_LO);

		return true;
	}


    int PowerCombinerParent::GetNElementsPerStrip()
    {
    	return fnElementsPerStrip;
    }

    void PowerCombinerParent::SetNElementsPerStrip( int aNumberOfElements )
    {
    	fnElementsPerStrip = aNumberOfElements;
    }

    double PowerCombinerParent::GetJunctionLoss()
    {
    	return fjunctionLoss;
    }
    void PowerCombinerParent::SetJunctionLoss( double aJunctionLoss )
    {
    	fjunctionLoss = aJunctionLoss;
    }
    double PowerCombinerParent::GetPatchLoss()
    {
    	return fpatchLoss;
    }
    void PowerCombinerParent::SetPatchLoss( double aPatchLoss )
    {
    	fpatchLoss = aPatchLoss;
    }
    double PowerCombinerParent::GetAmplifierLoss()
    {
    	return famplifierLoss;
    }
    void PowerCombinerParent::SetAmplifierLoss( double aAmplifierLoss )
    {
    	famplifierLoss = aAmplifierLoss;
    }
    double PowerCombinerParent::GetEndPatchLoss()
    {
    	return fendPatchLoss;
    }
    void PowerCombinerParent::SetEndPatchLoss( double aEndPatchLoss )
    {
    	fendPatchLoss = aEndPatchLoss;
    }
    double PowerCombinerParent::GetJunctionResistance()
    {
    	return fjunctionResistance;
    }
    void PowerCombinerParent::SetJunctionResistance( double aJunctionResistance )
    {
    	fjunctionResistance = aJunctionResistance;
    }
    double PowerCombinerParent::GetDampingFactor( int z_index )
    {
    	return fdampingFactors[z_index];
    }
    void PowerCombinerParent::SetDampingFactor (int z_index, double aDampingFactor )
    {
    	fdampingFactors[z_index] = aDampingFactor;
    }



} /* namespace locust */

