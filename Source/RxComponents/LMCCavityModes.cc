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
		fOrbitPhase( 0. ),
		fVoltagePhase( 0. )
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


        fRollingAvg.resize(GetNCavityModes());
        fCounter.resize(GetNCavityModes());
        for (int i = 0; i < GetNCavityModes(); i++)
        {
            fRollingAvg[i].resize(GetNCavityModes());
            fCounter[i].resize(GetNCavityModes());
            for (int j = 0; j < GetNCavityModes(); j++)
            {
            	fRollingAvg[i][j].resize(GetNCavityModes());
            	fCounter[i][j].resize(GetNCavityModes());
            }
        }

    	return true;
    }

	bool CavityModes::AddOneModeToCavityProbe(Signal* aSignal, std::vector<double> particleXP, double excitationAmplitude, double EFieldAtProbe, std::vector<double> cavityDopplerFrequency, double dt, double phi_LO, double totalScalingFactor, unsigned sampleIndex, bool initParticle)
	{
		double dopplerFrequency = cavityDopplerFrequency[0];
        SetVoltagePhase( GetVoltagePhase() + dopplerFrequency * dt ) ;
        double voltageValue = excitationAmplitude * EFieldAtProbe;
	//std::cout << "Amp, FieldAtProbe, Product: " << excitationAmplitude << " " << EFieldAtProbe << " " << voltageValue << std::endl;
        voltageValue *= cos(GetVoltagePhase());
	
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


			if ((l==GetNCavityModes()-1)&&(m==GetNCavityModes()-1)&&(n==GetNCavityModes()-1))
			{
				double totalEnergy = 0.;
				for (int iL=0; iL<GetNCavityModes(); iL++)
				{
					for (int iM=0; iM<GetNCavityModes(); iM++)
					{
						for (int iN=0; iN<GetNCavityModes(); iN++)
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

    double CavityModes::GetVoltagePhase()
    {
    	return fVoltagePhase;
    }
    void CavityModes::SetVoltagePhase ( double aPhase )
    {
        fVoltagePhase = aPhase;
    }


} /* namespace locust */
