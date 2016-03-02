/*
 * LMCBasebandTrack.cc
 *
 *  Created on: Jan 14 2015
 *      Author: plslocum after nsoblath
 */

#include "LMCBasebandTrackGenerator.hh"

#include "LMCLogger.hh"

using std::string;

namespace locust
{
    LMCLOGGER( lmclog, "BasebandTrackGenerator" );

    MT_REGISTER_GENERATOR(BasebandTrackGenerator, "baseband-track");

    BasebandTrackGenerator::BasebandTrackGenerator( const std::string& aName ) :
            Generator( aName ),
            fDoGenerateFunc( &BasebandTrackGenerator::DoGenerateFreq ),
            fElectronEnergy( 30.0 ), // keV
            fTotalLOFreqs( 26. ) // GHz
    {
        fRequiredSignalState = Signal::kFreq;
    }

    BasebandTrackGenerator::~BasebandTrackGenerator()
    {
    }

    bool BasebandTrackGenerator::Configure( const ParamNode* aParam )
    {
        if( aParam == NULL) return true;

        SetElectronEnergy( aParam->GetValue< double >( "electron-energy", fElectronEnergy ) );
        SetTotalLOFreqs( aParam->GetValue< double >( "total-lo-freqs", fTotalLOFreqs ) );


        if( aParam->Has( "domain" ) )
        {
            string domain = aParam->GetValue( "domain" );
            if( domain == "time" )
            {
                SetDomain( Signal::kTime );
                LMCDEBUG( lmclog, "Domain is equal to time.");
            }
            else if( domain == "freq" )
            {
                SetDomain( Signal::kFreq );
            }
            else
            {
                LMCERROR( lmclog, "Unable to use domain requested: <" << domain << ">" );
                return false;
            }
        }

        return true;
    }

    void BasebandTrackGenerator::Accept( GeneratorVisitor* aVisitor ) const
    {
        aVisitor->Visit( this );
        return;
    }

    double BasebandTrackGenerator::GetElectronEnergy() const
    {
        return fElectronEnergy;
    }

    void BasebandTrackGenerator::SetElectronEnergy( double aElectronEnergy )
    {
        fElectronEnergy = aElectronEnergy;
        return;
    }

    double BasebandTrackGenerator::GetTotalLOFreqs() const
    {
        return fTotalLOFreqs;
    }

    void BasebandTrackGenerator::SetTotalLOFreqs( double aTotalLOFreqs )
    {
        fTotalLOFreqs = aTotalLOFreqs;
        return;
    }

 

    Signal::State BasebandTrackGenerator::GetDomain() const
    {
        return fRequiredSignalState;
    }

    void BasebandTrackGenerator::SetDomain( Signal::State aDomain )
    {
        if( aDomain == fRequiredSignalState ) return;
        fRequiredSignalState = aDomain;  // pls changed == to =.
        if( fRequiredSignalState == Signal::kTime )
        {
            fDoGenerateFunc = &BasebandTrackGenerator::DoGenerateTime;
        }
        else if( fRequiredSignalState == Signal::kFreq )
        {
            fDoGenerateFunc = &BasebandTrackGenerator::DoGenerateFreq;
        }
        else
        {
            LMCWARN( lmclog, "Unknown domain requested: " << aDomain );
        }
        return;
    }

    double BasebandTrackGenerator::CalculateLarmorPower( double gamma ) const  
    {
    double power = 1./(4.*PI*8.85e-12)*(2./3.)*pow(1.602e-19,4.)/pow(9.11e-31,2.)/3.e8*
      pow(1.,2.)*(gamma*gamma-1.)*pow(sin(90.*PI/180.),2.)*(1./1.6027e-16);  // keV/s, B = 1 T, pitch angle = 90 deg.
    return power;
    }

    double BasebandTrackGenerator::CalculateGamma( double KineticEnergy ) const
    {
    double Gamma = 1. + KineticEnergy/511.;  // 511. keV
    return Gamma;
    }

    double BasebandTrackGenerator::CalculateCyclotronFrequency( double Gamma ) const
    {
    double Frequency = 1.602e-19*1./(2.*PI*Gamma*9.11e-31);  // Hz, B = 1 T
    return Frequency;
    }

    double BasebandTrackGenerator::CalculateBasebandFrequency( double CyclotronFrequency ) const
    {
    double BasebandFrequency = (CyclotronFrequency - fTotalLOFreqs*1.e9);  // Hz
    return BasebandFrequency;
    }

    bool BasebandTrackGenerator::DoGenerate( Signal* aSignal ) const
    {
        return (this->*fDoGenerateFunc)( aSignal );
    }

    bool BasebandTrackGenerator::DoGenerateTime( Signal* aSignal ) const
    {

        RunLengthCalculator *RunLengthCalculator1 = new RunLengthCalculator;
        double time = 0.;
        double TimeDependentEnergy = fElectronEnergy;  // keV, starting value.
        double LarmorPower = 0.;
        double TimeDependentAmplitude = 0.0; 
        double BBFreq = 0.;  // baseband frequency
        double ElectronStartTime = 0.001; // seconds
        double ElectronDuration = 0.001; // seconds

        for( unsigned index = 0; index < aSignal->TimeSize(); ++index )
        {
            time = (double)index/(RunLengthCalculator1->GetAcquisitionRate()*1.e6); // seconds
            if (time > ElectronStartTime && time< ElectronStartTime+ElectronDuration)  // typical track length in seconds.
              {
              LarmorPower = this->CalculateLarmorPower(this->CalculateGamma(TimeDependentEnergy));  // keV/s
              // Now lose LarmorPower*dt in keV          
              TimeDependentEnergy -= LarmorPower*RunLengthCalculator1->GetBinWidth(); 
              TimeDependentAmplitude = pow(LarmorPower*1.602e-16, 0.5); // keV/s * 1.6e-19 J/eV * 1.e3 eV/keV = W.
              BBFreq = this->CalculateBasebandFrequency(this->CalculateCyclotronFrequency(this->CalculateGamma(TimeDependentEnergy)));
//              printf("bbfreq is %g\n", BBFreq);
//              printf("cyclotron freq is %g\n", this->CalculateCyclotronFrequency(this->CalculateGamma(TimeDependentEnergy))); getchar();
//              printf("amplitude is %g\n", TimeDependentAmplitude*TimeDependentAmplitude); getchar();
              aSignal->SignalTime( index ) += TimeDependentAmplitude*cos(2.*PI*BBFreq*(time-ElectronStartTime));
              }
        }
        delete RunLengthCalculator1;
        return true;
    }

    bool BasebandTrackGenerator::DoGenerateFreq( Signal* aSignal ) const
    {
// This is not finished.
        RunLengthCalculator *RunLengthCalculator1 = new RunLengthCalculator;
        for( unsigned index = 0; index < aSignal->FreqSize(); ++index )
        {
            aSignal->SignalFreq( index )[0] += 0.1*cos(2.*3.1415926*4000.*(double)index/(RunLengthCalculator1->GetAcquisitionRate()*1.e6)); 
            aSignal->SignalFreq( index )[1] += 0.1*sin(2.*3.1415926*4000.*(double)index/(RunLengthCalculator1->GetAcquisitionRate()*1.e6)); 
        }
        delete RunLengthCalculator1;
        return true;
    }

} /* namespace locust */
