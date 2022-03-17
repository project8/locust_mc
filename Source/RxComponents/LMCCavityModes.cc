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

    CavityModes::CavityModes():
		fnCavityProbes( 0 ),
		fNCavityModes( 0 ),
		fCavityProbeInductance( 1.0 ),
		fCavityProbeZ( 0. ),
		fCavityProbeTheta( 0. )
    {
    }

    CavityModes::~CavityModes()
    {
    }


    bool CavityModes::Configure( const scarab::param_node& aParam )
    {

    	if( !PowerCombiner::Configure(aParam))
    	{
    		LERROR(lmclog,"Error configuring PowerCombiner class from CavityModes subclass");
    		return false;
    	}

    	if ( aParam.has( "cavity-probe-inductance" ) )
    	{
    		SetCavityProbeInductance(aParam["cavity-probe-inductance"]().as_double());
    	}

        if( aParam.has( "n-modes" ) )
        {
            fNCavityModes = aParam["n-modes"]().as_int();
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

    int CavityModes::GetNCavityProbes()
    {
        return fnCavityProbes;
    }

    void CavityModes::SetNCavityProbes( int aNumberOfProbes )
    {
     	fnCavityProbes = aNumberOfProbes;
    }

	bool CavityModes::AddOneModeToCavityProbe(Signal* aSignal, double excitationAmplitude, double dopplerFrequency, double dt, double phi_LO, double totalScalingFactor, unsigned sampleIndex)
	{

		double voltageValue = excitationAmplitude;

		aSignal->LongSignalTimeComplex()[sampleIndex][0] += 2. * voltageValue * totalScalingFactor * sin(phi_LO);
		aSignal->LongSignalTimeComplex()[sampleIndex][1] += 2. * voltageValue * totalScalingFactor * cos(phi_LO);


		if ( GetVoltageCheck() && (sampleIndex%100 < 1) )
			LPROG( lmclog, "Voltage " << sampleIndex << " is <" << aSignal->LongSignalTimeComplex()[sampleIndex][1] << ">" );
		return true;
	}


	bool CavityModes::AddOneSampleToRollingAvg(int l, int m, int n, double excitationAmplitude, unsigned sampleIndex)
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

    double CavityModes::GetCavityProbeInductance()
    {
    	return fCavityProbeInductance;
    }
    void CavityModes::SetCavityProbeInductance( double anInductance )
    {
    	fCavityProbeInductance = anInductance;
    }

    bool CavityModes::SetCavityProbeLocations(int nCavityProbes, double cavityLength)
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



    std::vector<double> CavityModes::GetCavityProbeZ()
    {
    	return fCavityProbeZ;
    }
    void CavityModes::SetCavityProbeZ ( std::vector<double> aVector )
    {
    	fCavityProbeZ = aVector;
    }
    std::vector<double> CavityModes::GetCavityProbeTheta()
    {
    	return fCavityProbeTheta;
    }
    void CavityModes::SetCavityProbeTheta ( std::vector<double> aVector )
    {
    	fCavityProbeTheta = aVector;
    }




} /* namespace locust */
