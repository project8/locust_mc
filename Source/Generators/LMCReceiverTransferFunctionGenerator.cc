/*
 * LMCReceiverTransferFunctionGenerator.cc
 *
 *  Created on: Feb 4, 2014
 *      Author: plslocum after nsoblath
 */

#include "LMCReceiverTransferFunctionGenerator.hh"

#include "logger.hh"

using std::string;

namespace locust
{
    LOGGER( lmclog, "ReceiverTransferFunctionGenerator" );

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

    	// low pass filter.
        for( unsigned index = 0; index < aSignal->FreqSize(); ++index )
        {
        	if (index > aSignal->FreqSize()*0.95) // shifted?
//            if (index > aSignal->FreqSize()/2*0.95 && index < aSignal->FreqSize()/2*1.05) // unshifted?
//            if (index > aSignal->FreqSize()*0.95) // shifted and real and positive?


        	{
//        		printf("aSignal %d is %g\n", index, aSignal->SignalFreq()[index][0]);
              aSignal->SignalFreq()[index][0] = 0.;
              aSignal->SignalFreq()[index][1] = 0.;
//      		printf("after lpf aSignal %d is %g\n\n", index, aSignal->SignalFreq()[index][0]);
//      		getchar();
        	}
        	else
        	{

            aSignal->SignalFreq()[index][0] /= (double)aSignal->FreqSize();
            aSignal->SignalFreq()[index][1] /= (double)aSignal->FreqSize();

        	}

        }
        return true;
    }

} /* namespace locust */
