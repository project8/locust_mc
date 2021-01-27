/*
 * LMC[name]Generator.cc
 *
 *  Created on: Mar 12, 2014
 *      Author: nsoblath
 */

#include "LMC[name]Generator.hh"

#include "logger.hh"

namespace locust
{
    LOGGER( lmclog, "[name]Generator" );

    MT_REGISTER_GENERATOR([name]Generator, "config-name");

    [name]Generator::[name]Generator( const std::string& aName ) :
//        fRF_frequency( 10.e6 ),
        Generator( aName ),
        fDoGenerateFunc( &[name]Generator::DoGenerateTime )
    {
        fRequiredSignalState = Signal::k[domain];
    }

    [name]Generator::~[name]Generator()
    {
    }

    bool [name]Generator::Configure( const scarab::param_node& aParam )
    {
/*
        if( aParam.has( "rf-frequency" ) )
        {
            fRF_frequency = aParam.get_value< double >( "rf-frequency", fRF_frequency );
        }
*/
        return true;
    }

    void [name]Generator::Accept( GeneratorVisitor* aVisitor ) const
    {
        aVisitor->Visit( this );
        return;
    }

    Signal::State [name]Generator::GetDomain() const
    {
        return fRequiredSignalState;
    }

    void [name]Generator::SetDomain( Signal::State aDomain )
    {
        if( aDomain == fRequiredSignalState ) return;
        fRequiredSignalState = aDomain;
        if( fRequiredSignalState == Signal::kTime )
        {
            fDoGenerateFunc = &[name]Generator::DoGenerateTime;
        }
        else if( fRequiredSignalState == Signal::kFreq )
        {
            fDoGenerateFunc = &[name]Generator::DoGenerateFreq;
        }
        else
        {
            LWARN( lmclog, "Unknown domain requested: " << aDomain );
        }
        return;
    }


    bool [name]Generator::DoGenerate( Signal* aSignal )
    {
        return (this->*fDoGenerateFunc)( aSignal );
    }

    bool [name]Generator::DoGenerateTime( Signal* aSignal )
    {
/*
        double RF_freq = 50.e6; // Hz
  	    double amplitude = 1.e-5; // volts
       
        for( unsigned index = 0; index < aSignal->TimeSize(); ++index )
        {
            int nRFHalfCycles = int(index*(1./(fAcquisitionRate*1.e6)*(2.0*RF_freq)));
//            int nRFHalfCycles = int(index*(1./(fAcquisitionRate*1.e6)*(2.0*fRF_frequency)));
            double voltage = amplitude; // square wave is high
	        if ( nRFHalfCycles%2==1 ) voltage = -amplitude; // square wave is low
            aSignal->SignalTimeComplex()[index][0] +=  voltage; // in-phase (I)
//            aSignal->SignalTimeComplex()[index][1] += [...]; // quadrature (Q) = 0 for real signals.
        }
*/
        return true;
    }

    bool [name]Generator::DoGenerateFreq( Signal* aSignal )
    {
        return true;
    }


} /* namespace locust */
