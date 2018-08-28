/*
 * LMCFakeTrackGenerator.cc
 *
 *  Created on: Aug 8 2018
 *      Author: plslocum
 */

#include <cmath>
#include "LMCFakeTrackGenerator.hh"
#include "LMCDigitizer.hh"
#include "logger.hh"
#include "LMCConst.hh"


using std::string;

namespace locust
{
    LOGGER( lmclog, "FakeTrackGenerator" );

    MT_REGISTER_GENERATOR(FakeTrackGenerator, "fake-track");

    FakeTrackGenerator::FakeTrackGenerator( const std::string& aName ) :
            Generator( aName ),
            fDoGenerateFunc( &FakeTrackGenerator::DoGenerateTime ),
            fSignalPower( 0. ),
            fStartFrequency( 0. ),
            fStartVPhase( 0. ),
            fSlope( 0. ),
            fStartTime( 0. ),
            fEndTime( 0. ),
            fLO_frequency( 0. )
    {
        fRequiredSignalState = Signal::kTime;
    }

    FakeTrackGenerator::~FakeTrackGenerator()
    {
    }


  bool FakeTrackGenerator::Configure( const scarab::param_node* aParam )
  {
    if( aParam == NULL) return true;

    if( aParam->has( "signal-power" ) )
    SetSignalPower( aParam->get_value< double >( "signal-power", fSignalPower ) );

    if( aParam->has( "start-frequency" ) )
    SetStartFrequency( aParam->get_value< double >( "start-frequency", fStartFrequency ) );

    if( aParam->has( "start-vphase" ) )
    SetStartVPhase( aParam->get_value< double >( "start-vphase", fStartVPhase ) );

    if( aParam->has( "slope" ) )
    SetSlope( aParam->get_value< double >( "slope", fSlope ) );

    if( aParam->has( "start-time" ) )
    SetStartTime( aParam->get_value< double >( "start-time", fStartTime ) );

    if( aParam->has( "end-time" ) )
    SetEndTime( aParam->get_value< double >( "end-time", fEndTime ) );

    if( aParam->has( "lo-frequency" ) )
    SetFrequency( aParam->get_value< double >( "lo-frequency", fLO_frequency ) );



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


    void FakeTrackGenerator::Accept( GeneratorVisitor* aVisitor ) const
    {
        aVisitor->Visit( this );
        return;
    }

    double FakeTrackGenerator::GetSignalPower() const
    {
        return fSignalPower;
    }

    void FakeTrackGenerator::SetSignalPower( double aPower )
    {
        fSignalPower = aPower;
        return;
    }

    double FakeTrackGenerator::GetStartFrequency() const
    {
        return fStartFrequency;
    }

    void FakeTrackGenerator::SetStartFrequency( double aFrequency )
    {
        fStartFrequency = aFrequency;
        return;
    }

    double FakeTrackGenerator::GetStartVPhase() const
    {
        return fStartVPhase;
    }

    void FakeTrackGenerator::SetStartVPhase( double aPhase )
    {
        fStartVPhase = aPhase;
        return;
    }

    double FakeTrackGenerator::GetSlope() const
    {
        return fSlope;
    }

    void FakeTrackGenerator::SetSlope( double aSlope )
    {
        fSlope = aSlope;
        return;
    }

    double FakeTrackGenerator::GetStartTime() const
    {
        return fStartTime;
    }

    void FakeTrackGenerator::SetStartTime( double aTime )
    {
        fStartTime = aTime;
        return;
    }

    double FakeTrackGenerator::GetEndTime() const
    {
        return fEndTime;
    }

    void FakeTrackGenerator::SetEndTime( double aTime )
    {
        fEndTime = aTime;
        return;
    }

    double FakeTrackGenerator::GetFrequency() const
    {
        return fLO_frequency;
    }

    void FakeTrackGenerator::SetFrequency( double aFrequency )
    {
        fLO_frequency = aFrequency;
        return;
    }


 

    Signal::State FakeTrackGenerator::GetDomain() const
    {
        return fRequiredSignalState;
    }

    void FakeTrackGenerator::SetDomain( Signal::State aDomain )
    {
        if( aDomain == fRequiredSignalState ) return;
        fRequiredSignalState = aDomain;  // pls changed == to =.
        if( fRequiredSignalState == Signal::kTime )
        {
            fDoGenerateFunc = &FakeTrackGenerator::DoGenerateTime;
        }
        else if( fRequiredSignalState == Signal::kFreq )
        {
            fDoGenerateFunc = &FakeTrackGenerator::DoGenerateFreq;
        }
        else
        {
            LWARN( lmclog, "Unknown domain requested: " << aDomain );
        }
        return;
    }


    bool FakeTrackGenerator::DoGenerate( Signal* aSignal )
    {
        return (this->*fDoGenerateFunc)( aSignal );
    }

    bool FakeTrackGenerator::DoGenerateTime( Signal* aSignal )
    {

        RunLengthCalculator *RunLengthCalculator1 = new RunLengthCalculator;
        
        const unsigned nchannels = fNChannels;
        double LO_phase = 0.;
        double dt = 1./aSignal->DecimationFactor()/(RunLengthCalculator1->GetAcquisitionRate()*1.e6);

        for (unsigned ch = 0; ch < nchannels; ++ch)
        {
        double voltage_phase = fStartVPhase;
        double track_frequency = fStartFrequency;

        for( unsigned index = 0; index < aSignal->TimeSize()*aSignal->DecimationFactor(); ++index )

        {
        double time = (double)index/aSignal->DecimationFactor()/(RunLengthCalculator1->GetAcquisitionRate()*1.e6);

    	LO_phase += 2.*LMCConst::Pi()*fLO_frequency*dt;

        if ((time > fStartTime) && (time < fEndTime))
        {
            track_frequency += fSlope*1.e6/1.e-3*dt;
            voltage_phase += 2.*LMCConst::Pi()*track_frequency*(dt);


            aSignal->LongSignalTimeComplex()[ch*aSignal->TimeSize()*aSignal->DecimationFactor() + index][0] += sqrt(50.)*sqrt(fSignalPower)*cos(voltage_phase-LO_phase);
            aSignal->LongSignalTimeComplex()[ch*aSignal->TimeSize()*aSignal->DecimationFactor() + index][1] += sqrt(50.)*sqrt(fSignalPower)*cos(-LMCConst::Pi()/2. + voltage_phase-LO_phase);
        }

        }
        }
        delete RunLengthCalculator1;
        return true;
    }

    bool FakeTrackGenerator::DoGenerateFreq( Signal* aSignal )
    {
        return true;
    }

} /* namespace locust */
