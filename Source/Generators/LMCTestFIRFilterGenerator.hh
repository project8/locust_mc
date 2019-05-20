/*
 * LMCTestFIRFilterGenerator.hh
 *
 *  Created on: May 3 2019
 *      Author: plslocum
 */

#ifndef LMCTESTFIRFILTERGENERATOR_HH_
#define LMCTESTFIRFILTERGENERATOR_HH_

#include "LMCGenerator.hh"
#include "LMCRunLengthCalculator.hh"
#include "LMCFieldBuffer.hh"


namespace scarab
{
  class param_node;
}

namespace locust
{
  class Digitizer;

    /*!
     @class TestFIRFilterGenerator
     @author P. L. Slocum

     @brief Add Sine Wave to the signal.

     @details
     Can operate in time or frequency space

     Configuration name: "test-firfilter"

     Available configuration options:
     - "frequency": double -- Frequency of the sine wave.
     - "amplitude": double -- Amplitude of the sine wave.
     - "domain": string -- Determines whether the sinusoidal test signal is generated in the time 
            or frequency domain
    
     Available options: "time" and "freq" [default]

    */
    class TestFIRFilterGenerator : public Generator
    {
        public:
            TestFIRFilterGenerator( const std::string& aName = "test-firfilter" );
            virtual ~TestFIRFilterGenerator();

            bool Configure( const scarab::param_node* aNode );
            bool Configure2( const Digitizer* aDig );

            void Accept( GeneratorVisitor* aVisitor ) const;

            double GetRFFrequency() const;
            void SetRFFrequency( double aFrequency );

            double GetLOFrequency() const;
            void SetLOFrequency( double aFrequency );


            double GetAmplitude() const;
            void SetAmplitude( double aAmplitude );

            Signal::State GetDomain() const;
            void SetDomain( Signal::State aDomain );


        private:
            bool DoGenerate( Signal* aSignal );

            bool DoGenerateTime( Signal* aSignal );
            bool DoGenerateFreq( Signal* aSignal );

            bool (TestFIRFilterGenerator::*fDoGenerateFunc)( Signal* aSignal );
            double* GetFIRFilter(int nskips);
            int GetNFilterBins(double* filterarray);
            double GetFIRSample(double* filterarray, int nfilterbins, double dtfilter, unsigned channel);

            void InitializeBuffers(unsigned filterbuffersize, unsigned fieldbuffersize);
            void FillBuffers(double FieldAmplitude, double FieldPhase, double LOPhase, unsigned index, unsigned channel);
            void PopBuffers(unsigned channel);

            double* filterarray;
            double fRF_frequency;
            double fLO_frequency;
            double fAmplitude;
            double fFilter_resolution;
            std::string gfilter_filename;

            std::vector<std::deque<double>> EFieldBuffer;
            std::vector<std::deque<double>> EPhaseBuffer;
            std::vector<std::deque<double>> EAmplitudeBuffer;
            std::vector<std::deque<double>> EFrequencyBuffer;
            std::vector<std::deque<double>> LOPhaseBuffer;
            std::vector<std::deque<double>> IndexBuffer;
            std::vector<std::deque<double>> PatchFIRBuffer;
    };

} /* namespace locust */

#endif /* LMCTestFIRFilterGENERATOR_HH_ */

