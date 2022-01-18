/*
 * LMCPowerCombiner.cc
 *
 *  Created on: Feb 25, 2020
 *      Author: pslocum
 */

#include "LMCPowerCombiner.hh"
#include <iostream>
#include "logger.hh"

namespace locust
{

	LOGGER( lmclog, "PowerCombiner" );


    PowerCombiner::PowerCombiner():
			fnElementsPerStrip( 0 ),
			fdampingFactors( 0 ),
            fjunctionLoss( 1.0 ),
            fpatchLoss( 0.6 ),
            famplifierLoss( 0.66 ),
            fendPatchLoss( 1.0 ),
            fjunctionResistance( 0.3 ),
            fnCavityProbes( 0 ),
            fCavityProbeInductance( 1.0 ),
            fCavityProbeZ( 0. ),
            fCavityProbeTheta( 0. ),
			fvoltageCheck( false ),
			fRollingAvg( 0. )

    {}
    PowerCombiner::~PowerCombiner() {}


    bool PowerCombiner::Configure( const scarab::param_node& aParam )
    {

    	if ( aParam.has( "voltage-check" ) )
    	{
    		fvoltageCheck = aParam["voltage-check"]().as_bool();
    	}

        if( aParam.has( "nelements-per-strip" ) )
        {
            fnElementsPerStrip = aParam["nelements-per-strip"]().as_int();
            if (fnElementsPerStrip < 1)
            {
        		LERROR(lmclog,"PowerCombiner expects >= 1 elements per strip.");
            	return false;
            }
            fdampingFactors.resize( fnElementsPerStrip );
        }

        return true;

    }


    void PowerCombiner::SayHello()
    {
    	printf("powercombiner says hello\n"); getchar();
    }

    bool PowerCombiner::IsSinglePatch()
    {
    	return false;
    }

    Receiver* PowerCombiner::ChooseElement()
    {
    	PatchAntenna* aPatch = new PatchAntenna;
    	return aPatch;
    }


	bool PowerCombiner::AddOneVoltageToStripSum(Signal* aSignal, double VoltageFIRSample, double phi_LO, unsigned z_index, unsigned sampleIndex)
	{

		VoltageFIRSample *= GetDampingFactor(z_index);
		aSignal->LongSignalTimeComplex()[sampleIndex][0] += 2.*VoltageFIRSample * sin(phi_LO);
		aSignal->LongSignalTimeComplex()[sampleIndex][1] += 2.*VoltageFIRSample * cos(phi_LO);

		if ( (fvoltageCheck==true) && (sampleIndex%100 < 1) )
			LPROG( lmclog, "Voltage " << z_index << "  " << sampleIndex << " is <" << aSignal->LongSignalTimeComplex()[sampleIndex][1] << ">" );
		return true;
	}


	bool PowerCombiner::AddOneModeToCavityProbe(Signal* aSignal, double VoltageFIRSample, double phi_LO, double totalScalingFactor, double cavityProbeImpedance, unsigned sampleIndex)
	{

		aSignal->LongSignalTimeComplex()[sampleIndex][0] += 2. * VoltageFIRSample * totalScalingFactor * cavityProbeImpedance * sin(phi_LO);
		aSignal->LongSignalTimeComplex()[sampleIndex][1] += 2. * VoltageFIRSample * totalScalingFactor * cavityProbeImpedance * cos(phi_LO);

//		printf("signal is %g\n", aSignal->LongSignalTimeComplex()[sampleIndex][0]); getchar();

		if ( (fvoltageCheck==true) && (sampleIndex%100 < 1) )
			LPROG( lmclog, "Voltage " << sampleIndex << " is <" << aSignal->LongSignalTimeComplex()[sampleIndex][1] << ">" );
		return true;
	}

	void PowerCombiner::SetCounter( int aValue )
	{
		fCounter = aValue;
	}

	bool PowerCombiner::AddOneSampleToRollingAvg(double VoltageFIRSample, double phi_LO, double totalScalingFactor, double cavityProbeImpedance, unsigned sampleIndex)
	{
    	char buffer[60];
		double vI = 2. * VoltageFIRSample * totalScalingFactor * cavityProbeImpedance * sin(phi_LO);
		double vQ = 2. * VoltageFIRSample * totalScalingFactor * cavityProbeImpedance * cos(phi_LO);

		fRollingAvg = ( fRollingAvg * fCounter + vI*vI + vQ*vQ ) / ( fCounter + 1 );
		int a = sprintf(buffer, "output/modeEnergies.txt");
		const char *fpname = buffer;
		FILE *fp = fopen(fpname, "a");

		if ( (fCounter == 1000) )
		{
			fprintf(fp, "fRollingAvg is %g\n", fRollingAvg);
		}
		fCounter += 1;
		fclose (fp);

		return true;
	}






    int PowerCombiner::GetNElementsPerStrip()
    {
    	return fnElementsPerStrip;
    }

    void PowerCombiner::SetNElementsPerStrip( int aNumberOfElements )
    {
    	fnElementsPerStrip = aNumberOfElements;
    }

    int PowerCombiner::GetNCavityProbes()
    {
    	return fnCavityProbes;
    }

    void PowerCombiner::SetNCavityProbes( int aNumberOfProbes )
    {
    	fnCavityProbes = aNumberOfProbes;
    }

    double PowerCombiner::GetJunctionLoss()
    {
    	return fjunctionLoss;
    }
    void PowerCombiner::SetJunctionLoss( double aJunctionLoss )
    {
    	fjunctionLoss = aJunctionLoss;
    }
    double PowerCombiner::GetPatchLoss()
    {
    	return fpatchLoss;
    }
    void PowerCombiner::SetPatchLoss( double aPatchLoss )
    {
    	fpatchLoss = aPatchLoss;
    }
    double PowerCombiner::GetAmplifierLoss()
    {
    	return famplifierLoss;
    }
    void PowerCombiner::SetAmplifierLoss( double aAmplifierLoss )
    {
    	famplifierLoss = aAmplifierLoss;
    }
    double PowerCombiner::GetEndPatchLoss()
    {
    	return fendPatchLoss;
    }
    void PowerCombiner::SetEndPatchLoss( double aEndPatchLoss )
    {
    	fendPatchLoss = aEndPatchLoss;
    }
    double PowerCombiner::GetJunctionResistance()
    {
    	return fjunctionResistance;
    }
    void PowerCombiner::SetJunctionResistance( double aJunctionResistance )
    {
    	fjunctionResistance = aJunctionResistance;
    }
    double PowerCombiner::GetDampingFactor( int z_index )
    {
    	return fdampingFactors[z_index];
    }
    void PowerCombiner::SetDampingFactor (int z_index, double aDampingFactor )
    {
    	fdampingFactors[z_index] = aDampingFactor;
    }

    double PowerCombiner::GetCavityProbeInductance()
    {
    	return fCavityProbeInductance;
    }
    void PowerCombiner::SetCavityProbeInductance( double anInductance )
    {
    	fCavityProbeInductance = anInductance;
    }

    bool PowerCombiner::SetCavityProbeLocations(int nCavityProbes, double cavityLength)
    {

    	SetNCavityProbes(nCavityProbes);
    	std::vector<double> probeZ;
    	probeZ.resize(nCavityProbes);

    	std::vector<double> probeTheta;
    	probeTheta.resize(nCavityProbes);

    	double probeSpacing = cavityLength / ((double)nCavityProbes + 1.);

		for (unsigned index=0; index<probeZ.size(); index++)
		{
			probeZ[index] = -cavityLength/2. + (index+1)*probeSpacing;
			probeTheta[index] = 0.0;
		}

    	SetCavityProbeZ(probeZ);
    	SetCavityProbeTheta(probeTheta);

    	return true;
    }



    std::vector<double> PowerCombiner::GetCavityProbeZ()
    {
    	return fCavityProbeZ;
    }
    void PowerCombiner::SetCavityProbeZ ( std::vector<double> aVector )
    {
    	fCavityProbeZ = aVector;
    }
    std::vector<double> PowerCombiner::GetCavityProbeTheta()
    {
    	return fCavityProbeTheta;
    }
    void PowerCombiner::SetCavityProbeTheta ( std::vector<double> aVector )
    {
    	fCavityProbeTheta = aVector;
    }



} /* namespace locust */

