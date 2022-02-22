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
			fNCavityModes( 0 ),
            fCavityProbeInductance( 1.0 ),
            fCavityProbeZ( 0. ),
            fCavityProbeTheta( 0. ),
			fvoltageCheck( false ),
			fRollingAvg( 0. ),
			fCounter( 0 ),
			fVoltagePhase( 0. )

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

        fRollingAvg.resize(fNCavityModes);
        fCounter.resize(fNCavityModes);
        for (int i = 0; i < fNCavityModes; i++)
        {
            fRollingAvg[i].resize(fNCavityModes);
            fCounter[i].resize(fNCavityModes);
            for (int j = 0; j < fNCavityModes; j++)
            {
            	fRollingAvg[i][j].resize(fNCavityModes);
            	fCounter[i][j].resize(fNCavityModes);
            }
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

    void PowerCombiner::SetNCavityModes( int aNumberOfModes )
    {
    	fNCavityModes = aNumberOfModes;
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


	bool PowerCombiner::AddOneModeToCavityProbe(Signal* aSignal, double excitationAmplitude, double dopplerFrequency, double dt, double phi_LO, double totalScalingFactor, double cavityProbeImpedance, unsigned sampleIndex)
	{

		fVoltagePhase += dopplerFrequency * dt;
		double voltageValue = excitationAmplitude * cos(fVoltagePhase);

		aSignal->LongSignalTimeComplex()[sampleIndex][0] += 2. * voltageValue * totalScalingFactor * cavityProbeImpedance * sin(phi_LO);
		aSignal->LongSignalTimeComplex()[sampleIndex][1] += 2. * voltageValue * totalScalingFactor * cavityProbeImpedance * cos(phi_LO);

		if ( (fvoltageCheck==true) && (sampleIndex%100 < 1) )
			LPROG( lmclog, "Voltage " << sampleIndex << " is <" << aSignal->LongSignalTimeComplex()[sampleIndex][1] << ">" );
		return true;
	}

	bool PowerCombiner::AddOneSampleToRollingAvg(int l, int m, int n, double excitationAmplitude, unsigned sampleIndex)
	{
    	char buffer[60];
		double amp = excitationAmplitude;  // Kass electron current * J\cdot E, with optional resonance if !fBypassTF.

		fRollingAvg[l][m][n] = ( fRollingAvg[l][m][n] * fCounter[l][m][n] + pow(amp,2.) ) / ( fCounter[l][m][n] + 1 );
		int a = sprintf(buffer, "output/modeEnergies.txt");
		const char *fpname = buffer;
		FILE *fp = fopen(fpname, "a");

		if ( (sampleIndex%1000 < 1) )
		{
			printf("Writing to file:  sampleIndex is %d, fCounter is %d\n",
					sampleIndex, fCounter[l][m][n]);

			fprintf(fp, "%d%d%d %g\n", l, m, n, fRollingAvg[l][m][n]);


			if ((l==fNCavityModes-1)&&(m==fNCavityModes-1)&&(n==fNCavityModes-1))
			{
				double totalEnergy = 0.;
				for (int iL=0; iL<fNCavityModes; iL++)
				{
					for (int iM=0; iM<fNCavityModes; iM++)
					{
						for (int iN=0; iN<fNCavityModes; iN++)
						{
							if (!isnan(fRollingAvg[iL][iM][iN]))
							{
								totalEnergy += fRollingAvg[iL][iM][iN];
							}
						}
					}
				}

				fprintf(fp, "\ntotal energy is %g\n\n\n", totalEnergy);

			}

		}

		fCounter[l][m][n] += 1;
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

