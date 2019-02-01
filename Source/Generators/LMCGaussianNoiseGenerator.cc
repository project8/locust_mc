/*
 * LMCGaussianNoiseGenerator.cc
 *
 *  Created on: Feb 4, 2014
 *      Author: nsoblath
 */

#include "LMCGaussianNoiseGenerator.hh"

#include "logger.hh"


using std::string;

namespace locust
{
    LOGGER( lmclog, "GaussianNoiseGenerator" );

    MT_REGISTER_GENERATOR(GaussianNoiseGenerator, "gaussian-noise");

    GaussianNoiseGenerator::GaussianNoiseGenerator( const std::string& aName ) :
        Generator( aName ),
        fDoGenerateFunc( &GaussianNoiseGenerator::DoGenerateFreq ),
        fMean( 0. ),
        fSigma( 1. ),
        fUniDist( 0., 360. ),
        fNormDist( fMean, fSigma )
    {
        fRequiredSignalState = Signal::kFreq;
    }

    GaussianNoiseGenerator::~GaussianNoiseGenerator()
    {
    }

    bool GaussianNoiseGenerator::Configure( const scarab::param_node* aParam )
    {
        if( aParam == NULL) return true;


        if( aParam->has( "noise-floor" ) )
        {
//            tSigma = sqrt( aParam->get_value< double >( "noise-floor" ) * fAcquisitionRate * 1.e6);  // sampling rate
            fSigma = sqrt( aParam->get_value< double >( "noise-floor" ));
//            printf("sigma is %g and fAcquisitionRate is %g\n", tSigma, fAcquisitionRate); getchar();
        }
        else
        {
            fSigma = aParam->get_value< double >( "sigma");
        }

  //      SetMeanAndSigma( aParam->get_value< double >( "mean", fMean ), tSigma );

        if( aParam->has( "domain" ) )
        {
            string domain = aParam->get_value( "domain" );
            if( domain == "time" )
            {
                SetDomain( Signal::kTime );
                LDEBUG( lmclog, "Domain is equal to time.");
            }
            else if( domain == "freq" )
            {
                SetDomain( Signal::kFreq );
            }
            else
            {
                LERROR( lmclog, "Unable to use domain requested: <" << domain << ">" );
                return false;
            }
        }

        return true;
    }

    void GaussianNoiseGenerator::Accept( GeneratorVisitor* aVisitor ) const
    {
        aVisitor->Visit( this );
        return;
    }

    double GaussianNoiseGenerator::GetMean() const
    {
        return fMean;
    }

    void GaussianNoiseGenerator::SetMean( double aMean )
    {
        SetMeanAndSigma( aMean, fSigma );
        return;
    }

    double GaussianNoiseGenerator::GetSigma() const
    {
        return fSigma;
    }

    void GaussianNoiseGenerator::SetSigma( double aSigma )
    {
        SetMeanAndSigma( fMean, aSigma );
        return;
    }

    void GaussianNoiseGenerator::SetMeanAndSigma( double aMean, double aSigma )
    {
        fNormDist = std::normal_distribution< double >( aMean, aSigma );
        fUniDist = std::uniform_real_distribution< double >(0.,360.);
        fMean = aMean;
        fSigma = aSigma;
        return;
    }

    Signal::State GaussianNoiseGenerator::GetDomain() const
    {
        return fRequiredSignalState;
    }

    void GaussianNoiseGenerator::SetDomain( Signal::State aDomain )
    {
        if( aDomain == fRequiredSignalState ) return;
        fRequiredSignalState = aDomain;  // pls changed == to =.
        if( fRequiredSignalState == Signal::kTime )
        {
            fDoGenerateFunc = &GaussianNoiseGenerator::DoGenerateTime;
        }
        else if( fRequiredSignalState == Signal::kFreq )
        {
            fDoGenerateFunc = &GaussianNoiseGenerator::DoGenerateFreq;
        }
        else
        {
            LWARN( lmclog, "Unknown domain requested: " << aDomain );
        }
        return;
    }


    bool GaussianNoiseGenerator::DoGenerate( Signal* aSignal )
    {
        return (this->*fDoGenerateFunc)( aSignal );
    }

    bool GaussianNoiseGenerator::DoGenerateTime( Signal* aSignal )
    {

    	SetMeanAndSigma( fMean, fSigma * sqrt(fAcquisitionRate * 1.e6) );

        double gain=1.;
        const unsigned nchannels = fNChannels;
        double phi = 0.;  // voltage phase
        double mag = 0.;  // voltage mag

        for (int ch=0; ch<nchannels; ch++)
        {
            for( unsigned index = 0; index < aSignal->TimeSize(); ++index )
            {
                phi = fUniDist( fRNG );
                mag = fNormDist( fRNG );
                aSignal->SignalTimeComplex()[ch*aSignal->TimeSize() + index][0] += gain*sqrt(50.)* mag * cos(phi*LMCConst::Pi()/180.);
                aSignal->SignalTimeComplex()[ch*aSignal->TimeSize() + index][1] += gain*sqrt(50.)* mag * sin(phi*LMCConst::Pi()/180.);
                //	    printf("noise signal is %g\n", aSignal->SignalTimeComplex()[ch*aSignal->TimeSize() + index][1]); getchar();
            }
        }

        return true;
    }

    bool GaussianNoiseGenerator::DoGenerateFreq( Signal* aSignal )
    {
        for( unsigned index = 0; index < aSignal->FreqSize(); ++index )
        {
            aSignal->SignalFreq()[index][0] += fNormDist( fRNG );
            aSignal->SignalFreq()[index][1] += fNormDist( fRNG );
        }
        return true;
    }

} /* namespace locust */
