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
		fVoltagePhase( 0 ),
		fp( 0 ),
		fInterface( KLInterfaceBootstrapper::get_instance()->GetInterface() )
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

        SetNCavityModes(fInterface->fField->GetNModes());

        SizeNChannels(GetNChannels());


        if( aParam.has( "norm-check" ) )
        {
            std::string sFileName = GetOutputPath()+"/ModeEnergies.root";
            std::string sTextFileName = GetOutputPath()+"/ModeEnergies.txt";
            fp = fopen(sTextFileName.c_str(), "w");
            fclose(fp);

#ifdef ROOT_FOUND
            fRootHistoWriter = RootHistoWriter::get_instance();
            fRootHistoWriter->SetFilename(sFileName);
            fRootHistoWriter->OpenFile("RECREATE");
            fRootHistoWriter->CloseFile();
#endif
        }

        fRollingAvg.resize(2);
        fCounter.resize(2);	

	for ( int bTE = 0; bTE<2; bTE++)
	{
        	fRollingAvg[bTE].resize(GetNCavityModes());
        	fCounter[bTE].resize(GetNCavityModes());
        	for (int i = 0; i < GetNCavityModes(); i++)
        	{
            		fRollingAvg[bTE][i].resize(GetNCavityModes());
            		fCounter[bTE][i].resize(GetNCavityModes());
            		for (int j = 0; j < GetNCavityModes(); j++)
            		{
            			fRollingAvg[bTE][i][j].resize(GetNCavityModes());
            			fCounter[bTE][i][j].resize(GetNCavityModes());
			}
            }
        }


    	return true;
    }

    bool CavityModes::SizeNChannels(int aNumberOfChannels)
    {
    	SetNChannels(aNumberOfChannels);

    	std::vector<std::vector<std::vector<std::vector<double>>>> tZeroVector;
    	fVoltagePhase.swap(tZeroVector);
    	fVoltagePhase.resize(aNumberOfChannels);

    	for (int n = 0; n < GetNChannels(); n++)
    	{
    		fVoltagePhase[n].resize(GetNCavityModes());
    		for (int i = 0; i < GetNCavityModes(); i++)
    		{
    			fVoltagePhase[n][i].resize(GetNCavityModes());
    			for (int j = 0; j < GetNCavityModes(); j++)
    			{
    				fVoltagePhase[n][i][j].resize(GetNCavityModes());
    			}
    		}
    	}

    	return true;
    }

	bool CavityModes::AddOneModeToCavityProbe(int l, int m, int n, Signal* aSignal, std::vector<double> particleXP, double excitationAmplitude, double EFieldAtProbe, std::vector<double> cavityDopplerFrequency, double dt, double phi_LO, double totalScalingFactor, unsigned sampleIndex, int channelIndex, bool initParticle)
	{
		double dopplerFrequency = cavityDopplerFrequency[0];  // Only one shift, unlike in waveguide.
		SetVoltagePhase( GetVoltagePhase(channelIndex, l, m, n) + dopplerFrequency * dt, channelIndex, l, m, n ) ;
		double voltageValue = excitationAmplitude * EFieldAtProbe;
		voltageValue *= cos(GetVoltagePhase(channelIndex, l, m, n));

		aSignal->LongSignalTimeComplex()[sampleIndex][0] += 2. * voltageValue * totalScalingFactor * sin(phi_LO);
		aSignal->LongSignalTimeComplex()[sampleIndex][1] += 2. * voltageValue * totalScalingFactor * cos(phi_LO);

		if ( GetVoltageCheck() && (sampleIndex%100 < 1) )
			LPROG( lmclog, "Voltage " << sampleIndex << " is <" << aSignal->LongSignalTimeComplex()[sampleIndex][1] << ">" );
		return true;
	}

	bool CavityModes::WriteRootHisto()
	{
#ifdef ROOT_FOUND
		int nModes = GetNCavityModes()*100 + GetNCavityModes()*10 + GetNCavityModes();
		fRootHistoWriter->OpenFile("UPDATE");
		TH1D* aHisto = new TH1D("ModeEnergies", "TEMode Energy Depositions; Mode l*100 + m*10 + n; log10(Energy[arb])", nModes, 0., nModes);
		int bTE = 1;
		for (int iBin=0; iBin<nModes; iBin++)
		{
			// initialize histo:
			aHisto->SetBinContent(iBin+1, -100.);
		}

		for (int iL=0; iL<GetNCavityModes(); iL++)
		{
			for (int iM=0; iM<GetNCavityModes(); iM++)
			{
				for (int iN=0; iN<GetNCavityModes(); iN++)
				{
					int binIndex = iL*100 + iM*10 + iN;
					if (fRollingAvg[bTE][iL][iM][iN] > 0.)
					{
						aHisto->SetBinContent(binIndex+1, log10(fRollingAvg[bTE][iL][iM][iN]));
					}
				}
			}
		}

		fRootHistoWriter->Write1DHisto(aHisto);
		fRootHistoWriter->CloseFile();

#endif
		return true;
	}


	bool CavityModes::AddOneSampleToRollingAvg(int bTE, int l, int m, int n, double excitationAmplitude, unsigned sampleIndex)
	{
	    char cBufferFileName[60];
	    int a = sprintf(cBufferFileName, "%s/ModeEnergies.txt", GetOutputPath().c_str());
	    const char *cFileName = cBufferFileName;
	    fp = fopen(cFileName, "a");

	    double amp = excitationAmplitude;

	    fRollingAvg[bTE][l][m][n] = ( fRollingAvg[bTE][l][m][n] * fCounter[bTE][l][m][n] + pow(amp,2.) ) / ( fCounter[bTE][l][m][n] + 1 );

	    if ( (sampleIndex%1000 < 1) && (sampleIndex < 20000) )
	    {

		    fprintf(fp, "%d %d%d%d %g\n", bTE, l, m, n, fRollingAvg[bTE][l][m][n]);

		    if ((l==GetNCavityModes()-1)&&(m==GetNCavityModes()-1)&&(n==GetNCavityModes()-1))
		    {
			    double totalEnergy = 0.;
			    for (int iL=0; iL<GetNCavityModes(); iL++)
			    {
				    for (int iM=0; iM<GetNCavityModes(); iM++)
				    {
					    for (int iN=0; iN<GetNCavityModes(); iN++)
					    {
						    if (!isnan(fRollingAvg[bTE][iL][iM][iN]))
						    {
							    totalEnergy += fRollingAvg[bTE][iL][iM][iN];
						    }
					    }
				    }
			    }

			    fprintf(fp, "\ntotal energy is %g\n\n\n", totalEnergy);

#ifdef ROOT_FOUND
	    		    WriteRootHisto();
#endif
	    		    LPROG( lmclog, "\n\n\nMode energies written to files output/modeEnergies.root and output/modeEnergies.txt: sampleIndex is " << sampleIndex);
	    		    LPROG( lmclog, "\n\n\nPress return to continue averaging and writing mode energies, or Cntrl-C to quit.");
	    		    getchar();
	    	}

	    }

		fCounter[bTE][l][m][n] += 1;
		fclose (fp);

	    return true;
	}

    double CavityModes::GetVoltagePhase(int aChannel, int l, int m, int n)
    {
    	return fVoltagePhase[aChannel][l][m][n];
    }

    void CavityModes::SetVoltagePhase ( double aPhase, int aChannel, int l, int m, int n )
    {
        fVoltagePhase[aChannel][l][m][n] = aPhase;
    }


} /* namespace locust */
