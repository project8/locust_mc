/*
 * LMCLocalOscillatorGenerator.cc
 *
 *  Created on: Apr 20, 2018
 *      Author: buzinsky
 */

#include "LMCLocalOscillatorGenerator.hh"

#include "LMCRunLengthCalculator.hh"
#include "LMCSimulationController.hh"
#include "LMCConst.hh"

#include "logger.hh"

using std::string;

namespace locust
{
    LOGGER( lmclog, "LocalOscillatorGenerator" );

    MT_REGISTER_GENERATOR(LocalOscillatorGenerator, "local-oscillator");

    LocalOscillatorGenerator::LocalOscillatorGenerator( const std::string& aName ) :
            Generator( aName ),
            fDoGenerateFunc( &LocalOscillatorGenerator::DoGenerateTime )
    {
        fRequiredSignalState = Signal::kTime;
    }

    LocalOscillatorGenerator::~LocalOscillatorGenerator()
    {
    }

    bool LocalOscillatorGenerator::Configure( const scarab::param_node* aParam )
    {
        if( aParam == NULL) return true;

        if( aParam->has( "LO-Frequency" ) )
        {
            fLOFrequency = aParam->get_value< double >( "LO-Frequency" );
        }

        return true;
    }

    void LocalOscillatorGenerator::Accept( GeneratorVisitor* aVisitor ) const
    {
        aVisitor->Visit( this );
        return;
    }


    bool LocalOscillatorGenerator::DoGenerate( Signal* aSignal )
    {
        return (this->*fDoGenerateFunc)( aSignal );
    }

    bool LocalOscillatorGenerator::DoGenerateTime( Signal* aSignal )
    {
    	SimulationController SimulationController1;
        const unsigned nChannels = SimulationController1.GetNChannels();

        RunLengthCalculator *simpleCalculator =  new RunLengthCalculator();
        const double acquisitionTimeStep = 1. / ( simpleCalculator->GetAcquisitionRate() * 1.e6); //Rate is in MHz

        const unsigned signalSize = aSignal->TimeSize();
        const unsigned decimationFactor = aSignal->DecimationFactor();
        double *downmixedSignal;


        for(unsigned channelIndex = 0; channelIndex < nChannels; ++channelIndex)
        {
            double localOscillatorPhase = 0.;

            for( unsigned index = 0; index < decimationFactor*signalSize; ++index )
            {
                localOscillatorPhase += 2. * LMCConst::Pi() * fLOFrequency * acquisitionTimeStep;
                fftw_complex LOMixingPhase = {cos(localOscillatorPhase), -sin(localOscillatorPhase) };
                    ComplexMultiplication( aSignal->LongSignalTimeComplex()[channelIndex*signalSize*decimationFactor + index] , LOMixingPhase);

                //Assign to aSignal
                aSignal->LongSignalTimeComplex()[channelIndex*signalSize*decimationFactor + index][0] = downmixedSignal[0];
                aSignal->LongSignalTimeComplex()[channelIndex*signalSize*decimationFactor + index][1] = downmixedSignal[1];

            }
        }

        //Cleanup
        delete []downmixedSignal;

        return true;
    }

    double* LocalOscillatorGenerator::ComplexMultiplication(fftw_complex a, fftw_complex b)
    {
        double *product = new double[2];
        product[0] = a[0] * b[0] - a[1] * b[1];
        product[1] = a[1] * b[0] + a[0] * b[1];
        return product;
    }

} /* namespace locust */
