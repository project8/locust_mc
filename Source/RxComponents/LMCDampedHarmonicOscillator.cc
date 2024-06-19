/*
 * LMCDampedHarmonicOscillator.cc
 *
 *  Created on: Jul 6, 2022
 *      Author: pslocum
 */

#include "LMCDampedHarmonicOscillator.hh"
#include "logger.hh"

namespace locust
{

	LOGGER( lmclog, "DampedHarmonicOscillator" );

    DampedHarmonicOscillator::DampedHarmonicOscillator():
    		fMaxNBins( 20000 ),
			fTimeResolutionDefault( 1.e-10 ),
			fCavityFrequencyDefault( 1.067e9 ),
			fCavityQDefault( 1000 ),
			fThresholdFactorDefault ( 0.25 ),
			fHannekePowerFactorDefault( 1. )

    {}
    DampedHarmonicOscillator::~DampedHarmonicOscillator() {}

    bool DampedHarmonicOscillator::Configure( const scarab::param_node& aParam )
    {


    	if( !AnalyticResponseFunction::Configure(aParam))
    	{
    		LERROR(lmclog,"Error configuring AnalyticResponseFunction class from DampedHarmonicOscillator subclass");
    	}

    	if( aParam.has( "n-modes" ) )
    	{
    		SetNModes( aParam["n-modes"]().as_int() );
    	}
    	if( aParam.has( "dho-max-nbins" ) )
    	{
    		fMaxNBins = aParam["dho-max-nbins"]().as_int();
    	}
    	if( aParam.has( "dho-time-resolution" ) )
    	{
    		SetDHOTimeResolution( aParam["dho-time-resolution"]().as_double() );
    	}
    	if( aParam.has( "dho-threshold-factor" ) )
    	{
    		SetDHOThresholdFactor( aParam["dho-threshold-factor"]().as_double() );
    	}
    	if( aParam.has( "dho-cavity-frequency" ) )
    	{
    		SetCavityFrequency( aParam["dho-cavity-frequency"]().as_double() );
    	}
    	if( aParam.has( "dho-cavity-Q" ) )
    	{
    		SetCavityQ( aParam["dho-cavity-Q"]().as_double() );
    	}

    	fTFReceiverHandler = new TFReceiverHandler;
    	if(!fTFReceiverHandler->Configure(aParam))
    	{
    		LERROR(lmclog,"Error configuring receiver FIRHandler class");
    		exit(-1);
    		return false;
    	}

    	if ( !Initialize( GetNModes() ) )
    	{
    	    LERROR(lmclog,"Error while initializing DampedHarmonicOscillator.");
    	}


    	return true;
    }

    bool DampedHarmonicOscillator::Initialize( int nModes )
    {

        fCavityFrequency.resize(2);
        fCavityQ.resize(2);
        fTimeResolution.resize(2);
        fThresholdFactor.resize(2);
        fBFactor.resize(2);
        fCavityDampingFactor.resize(2);
        fHannekePowerFactor.resize(2);
        fCavityOmega.resize(2);
        fCavityOmegaPrime.resize(2);
        fHannekePowerFactor.resize(2);

        for(int bTE=0; bTE<2; bTE++)
        {
            fCavityFrequency[bTE].resize(nModes);
            fCavityQ[bTE].resize(nModes);
            fTimeResolution[bTE].resize(nModes);
            fThresholdFactor[bTE].resize(nModes);
            fBFactor[bTE].resize(nModes);
            fCavityDampingFactor[bTE].resize(nModes);
            fHannekePowerFactor[bTE].resize(nModes);
            fCavityOmega[bTE].resize(2);
            fCavityOmegaPrime[bTE].resize(2);

            for(int l=0; l<nModes; l++)
            {
                fCavityFrequency[bTE][l].resize(nModes);
                fCavityQ[bTE][l].resize(nModes);
                fTimeResolution[bTE][l].resize(nModes);
                fThresholdFactor[bTE][l].resize(nModes);
                fBFactor[bTE][l].resize(nModes);
                fCavityDampingFactor[bTE][l].resize(nModes);
                fHannekePowerFactor[bTE][l].resize(nModes);
                fCavityOmega[bTE][l].resize(nModes);
                fCavityOmegaPrime[bTE][l].resize(nModes);

                for(int m=0; m<nModes; m++)
                {
                    fCavityFrequency[bTE][l][m].resize(nModes);
                    fCavityQ[bTE][l][m].resize(nModes);
                    fTimeResolution[bTE][l][m].resize(nModes);
                    fThresholdFactor[bTE][l][m].resize(nModes);
                    fBFactor[bTE][l][m].resize(nModes);
                    fCavityDampingFactor[bTE][l][m].resize(nModes);
                    fHannekePowerFactor[bTE][l][m].resize(nModes);
                    fCavityOmega[bTE][l][m].resize(nModes);
                    fCavityOmegaPrime[bTE][l][m].resize(nModes);

                    for(int n=0; n<nModes; n++)
                    {
                    	fCavityFrequency[bTE][l][m][n] = fCavityFrequencyDefault;
                    	fCavityQ[bTE][l][m][n] = fCavityQDefault;
                    	fTimeResolution[bTE][l][m][n] = fTimeResolutionDefault;
                    	fThresholdFactor[bTE][l][m][n] = fThresholdFactorDefault;
                        fHannekePowerFactor[bTE][l][m][n] = fHannekePowerFactorDefault; // TO-DO:  Make more configurable.
                    	fCavityOmega[bTE][l][m][n] = fCavityFrequency[bTE][l][m][n] * 2. * LMCConst::Pi();
                    	fCavityDampingFactor[bTE][l][m][n] = 0.;
                    	if ( fCavityQ[bTE][l][m][n] > 0. )
                    	{
                    		fCavityDampingFactor[bTE][l][m][n] = 1. / 2. / fCavityQ[bTE][l][m][n];
                    	}
                    	fBFactor[bTE][l][m][n] = fCavityDampingFactor[bTE][l][m][n] * fCavityOmega[bTE][l][m][n];
                    	fCavityOmegaPrime[bTE][l][m][n] = sqrt( fCavityOmega[bTE][l][m][n]*fCavityOmega[bTE][l][m][n] - fBFactor[bTE][l][m][n]*fBFactor[bTE][l][m][n] );
                    }
                }
            }
        }

    	if (!GenerateGreensFunction()) return false;

        return true;
    }


    void DampedHarmonicOscillator::SetCavityQ( double aQ )
    {
        fCavityQDefault = aQ;
    }
    double DampedHarmonicOscillator::GetCavityQ( int bTE, int l, int m, int n )
    {
    	return fCavityQ[bTE][l][m][n];
    }
    void DampedHarmonicOscillator::SetCavityFrequency( double aFrequency )
    {
        fCavityFrequencyDefault = aFrequency;
    }
    double DampedHarmonicOscillator::GetCavityFrequency( int bTE, int l, int m, int n )
    {
    	return fCavityFrequency[bTE][l][m][n];
    }
    void DampedHarmonicOscillator::SetDHOTimeResolution( double aTimeResolution )
    {
        fTimeResolutionDefault = aTimeResolution;
    }
    double DampedHarmonicOscillator::GetDHOTimeResolution( int bTE, int l, int m, int n )
    {
    	return fTimeResolution[bTE][l][m][n];
    }
    void DampedHarmonicOscillator::SetDHOThresholdFactor( double aThresholdFactor )
    {
        fThresholdFactorDefault = aThresholdFactor;
    }
    double DampedHarmonicOscillator::GetDHOThresholdFactor( int bTE, int l, int m, int n )
    {
    	return fThresholdFactor[bTE][l][m][n];
    }

    double DampedHarmonicOscillator::ExpDecayTerm( int bTE, int l, int m, int n, double t)
    {
    	double ExpDecayTerm = exp( -fBFactor[bTE][l][m][n] * t);
    	return ExpDecayTerm;
    }

    std::pair<double,double> DampedHarmonicOscillator::GreensFunction( int bTE, int l, int m, int n, double t)
    {

    	double tExpDecayTerm = ExpDecayTerm( bTE, l, m, n, t );
    	double GreensFunctionValueReal = fTimeResolution[bTE][l][m][n] * fHannekePowerFactor[bTE][l][m][n] * tExpDecayTerm * sin( fCavityOmegaPrime[bTE][l][m][n] * t);
    	double GreensFunctionValueImag = -1. * fTimeResolution[bTE][l][m][n] * fHannekePowerFactor[bTE][l][m][n] * tExpDecayTerm * cos( fCavityOmegaPrime[bTE][l][m][n] * t);

    	return std::make_pair(GreensFunctionValueReal,GreensFunctionValueImag);
    }


    double DampedHarmonicOscillator::NormFactor(int bTE, int l, int m, int n, double aDriveFrequency)
    {
        if (!fTFReceiverHandler->ConvertAnalyticGFtoFIR({{bTE,l,m,n}}, GetGFarray( {{bTE,l,m,n}} )))
        {
            LERROR(lmclog,"GF->FIR was not generated in DHO::NormFactor.");
            exit(-1);
            return false;
        }

        /* initialize time series */
        Signal* aSignal = new Signal();
        int N0 = GetGFarray( {{bTE, l, m, n}} )[0].size();
        aSignal->Initialize( N0 , 1 );
        double convolutionMag = 0.;

        for (unsigned i=0; i<1000; i++)  // time stamps
        {
            // populate time series and convolve it with the FIR filter
            PopulateCalibrationSignal(bTE, l, m, n, aSignal, N0, aDriveFrequency);

            std::pair<double,double> convolutionPair = fTFReceiverHandler->ConvolveWithComplexFIRFilterArray(bTE,l,m,n,SignalToDequeArray(bTE, l, m, n, aSignal));

            if (fabs(convolutionPair.first) > convolutionMag)
            {
                convolutionMag = convolutionPair.first;
            }
        } // i

        delete aSignal;

        return convolutionMag;

    }

	bool DampedHarmonicOscillator::PopulateCalibrationSignal(int bTE, int l, int m, int n, Signal* aSignal, int N0, double aDriveFrequency)
	{

        double voltage_phase = 0.;

        for( unsigned index = 0; index < N0; ++index )
        {
            voltage_phase = 2.*LMCConst::Pi()*aDriveFrequency*(double)index*fTimeResolution[bTE][l][m][n];

            aSignal->LongSignalTimeComplex()[index][0] = cos(voltage_phase);
            aSignal->LongSignalTimeComplex()[index][1] = cos(-LMCConst::Pi()/2. + voltage_phase);
        }

        return true;
	}

	std::deque<double> DampedHarmonicOscillator::SignalToDequeArray(int bTE, int l, int m, int n, Signal* aSignal)
	{
	    std::deque<double> incidentSignal;
	    for (unsigned i=0; i<fTFReceiverHandler->GetFilterSizeArray(bTE, l, m, n); i++)
	    {
	    	incidentSignal.push_back(aSignal->LongSignalTimeComplex()[i][0]);
	    }
	    return incidentSignal;
	}


    bool DampedHarmonicOscillator::GenerateGreensFunction()
    {
        int nModes = GetNModes();

        std::vector <std::vector< std::vector< std::vector< std::vector<std::pair<double,std::pair<double,double> > > > > > > tGFArray;
        tGFArray.resize(2); // TE/TM

        for( int bTE=0; bTE<2; bTE++)
        {
            tGFArray[bTE].resize(nModes);

            for(int l=0; l<nModes; l++)
            {
                tGFArray[bTE][l].resize(nModes);
                for(int m=0; m<nModes; m++)
                {
                    tGFArray[bTE][l][m].resize(nModes);
                    for(int n=0; n<nModes; n++)
                    {
                        for (unsigned i=0; i<fMaxNBins; i++)
                        {
                            double tValue = i * fTimeResolution[bTE][l][m][n];
                            if ( ExpDecayTerm( bTE, l, m, n, tValue) > fThresholdFactor[bTE][l][m][n] * ExpDecayTerm(bTE, l, m, n, 0.) )
                            {
                                tGFArray[bTE][l][m][n].push_back(std::make_pair(fTimeResolution[bTE][l][m][n],GreensFunction( bTE, l, m, n, tValue)));
                            }
                        }
                        std::reverse( tGFArray[bTE][l][m][n].begin(), tGFArray[bTE][l][m][n].end() );
                    }
                }
            }
        }

        SetGFarray(tGFArray ); // unnormalized

        for( int bTE=0; bTE<2; bTE++)
        {
            for(int l=0; l<nModes; l++)
            {
                for(int m=0; m<nModes; m++)
                {
                    for(int n=0; n<nModes; n++)
                    {
                        double aNormFactor = NormFactor( bTE, l, m, n, fCavityFrequency[bTE][l][m][n]);
                        for (unsigned i=0; i<tGFArray[bTE][l][m][n].size(); i++)
                        {
                            if (aNormFactor > 0.)
                            {
                                tGFArray[bTE][l][m][n][i].second.first /= aNormFactor;
                                tGFArray[bTE][l][m][n][i].second.second /= aNormFactor;
                            }
                        }
                    }
                }
            }
        }

        SetGFarray( tGFArray); // now normalized.

    	if ( tGFArray.size() < 1 ) return false;
    	else return true;

    }




} /* namespace locust */

