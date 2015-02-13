/*
 * LMCReceiverTransferFunctionGenerator.cc
 *
 *  Created on: Feb 4, 2014
 *      Author: plslocum after nsoblath
 */

#include "LMCReceiverTransferFunctionGenerator.hh"

#include "../Core/LMCLogger.hh"

using std::string;

namespace locust
{
    LMCLOGGER( lmclog, "ReceiverTransferFunctionGenerator" );

    MT_REGISTER_GENERATOR(ReceiverTransferFunctionGenerator, "receiver-transfer-function");

    ReceiverTransferFunctionGenerator::ReceiverTransferFunctionGenerator( const std::string& aName ) :
            Generator( aName ),
            fDoGenerateFunc( &ReceiverTransferFunctionGenerator::DoGenerateFreq ),
            fReceiverGain( 0. ) // dB
    {
        fRequiredSignalState = Signal::kFreq;
    }

    ReceiverTransferFunctionGenerator::~ReceiverTransferFunctionGenerator()
    {
    }

    bool ReceiverTransferFunctionGenerator::Configure( const ParamNode* aParam )
    {
        if( aParam == NULL) return true;
        SetReceiverGain( aParam->GetValue< double >( "receiver-gain", fReceiverGain ) );
        return true;
    }

    void ReceiverTransferFunctionGenerator::Accept( GeneratorVisitor* aVisitor ) const
    {
        aVisitor->Visit( this );
        return;
    }


    double ReceiverTransferFunctionGenerator::GetReceiverGain() const
    {
        return fReceiverGain;
    }

    void ReceiverTransferFunctionGenerator::SetReceiverGain( double aReceiverGain )
    {
        fReceiverGain = aReceiverGain;
        return;
    }


    bool ReceiverTransferFunctionGenerator::DoGenerate( Signal* aSignal ) const
    {
        return (this->*fDoGenerateFunc)( aSignal );
    }

    bool ReceiverTransferFunctionGenerator::DoGenerateTime( Signal* aSignal ) const
    {
        // Nothing happens to the signal if this generator operates in the time domain.
        return true;
    }



// This function needs to be checked.  It seems to be applying a factor 
//to the frequency spectrum even when the transfer function == 1.0 .
    bool ReceiverTransferFunctionGenerator::DoGenerateFreq( Signal* aSignal ) const
    {
        double VoltageGain = pow(10.,fReceiverGain/2./10.);  // divide by 2 if gain is given in power.
        double SqrtImpedance = pow(1., 0.5);  //  assume 1 ohm impedance for simple testing.

//        for( unsigned index = 0; index < aSignal->FreqSize(); ++index )
//        {

        // Get through antenna and apply voltage gain.  Assume perfect coupling.
//            aSignal->SignalFreq( index )[0] *= VoltageGain * SqrtImpedance;
//            aSignal->SignalFreq( index )[1] *= VoltageGain * SqrtImpedance;

        // Step function for testing.  Factor of 1.0
//            if (index<aSignal->FreqSize()/2.)
//            {  
//            aSignal->SignalFreq( index )[0] *= 1.0;
//            aSignal->SignalFreq( index )[1] *= 1.0;
//            }
//            else
//            {  
//            aSignal->SignalFreq( index )[0] *= 1.0;
//            aSignal->SignalFreq( index )[1] *= 1.0;
//            }
//        }
        return true;
    }

} /* namespace locust */
