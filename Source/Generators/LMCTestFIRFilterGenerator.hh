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
#include "LMCHilbertTransform.hh"


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

     @brief Compute FIR response to a sinusoidal E field, treating the E field as a real arbitrary signal.  Buffer E fields to constrain arrival times (roughly).  Extract mag and phase of incident E field with Hilbert transform and convolve it with FIR filter to generate each voltage.  

     @details
     Can operate in time or frequency space

     Configuration name: "test-firfilter"

     Available configuration options:
     - "rf-frequency": double -- Frequency of the incident sine wave.
     - "lo-frequency": double -- Frequency of the local oscillator.
     - "amplitude": double -- Amplitude of the incident sine wave.
     - "filter-filename": double -- path to FIR text file.
     - "filter-resolution": double -- time resolution of coefficients in filter-filename.
     - "buffer-size": double -- size of buffer to contain incident E field values.
     - "buffer-margin": double -- distance from beginning of buffer at which to extract Hilbert transform info.  extrapolate to edge across this margin.
     - "domain": string -- Determines whether the sinusoidal test signal is generated in the time 
            or frequency domain
    
     Available options: "time" and "freq" [default]

    */
    class TestFIRFilterGenerator : public Generator
    {
        public:
            TestFIRFilterGenerator( const std::string& aName = "test-firfilter" );
            virtual ~TestFIRFilterGenerator();

            bool Configure( const scarab::param_node& aNode );
            bool Configure2( const Digitizer* aDig );

            void Accept( GeneratorVisitor* aVisitor ) const;

            double GetRFFrequency() const;
            void SetRFFrequency( double aFrequency );

            double GetLOFrequency() const;
            void SetLOFrequency( double aFrequency );

            double GetAmplitude() const;
            void SetAmplitude( double aAmplitude );

            double GetBufferSize() const;
            void SetBufferSize( double aBufferSize );

            double GetBufferMargin() const;
            void SetBufferMargin( double aBufferMargin );

            Signal::State GetDomain() const;
            void SetDomain( Signal::State aDomain );


        private:
            bool DoGenerate( Signal* aSignal );

            bool DoGenerateTime( Signal* aSignal );
            bool DoGenerateFreq( Signal* aSignal );

            bool (TestFIRFilterGenerator::*fDoGenerateFunc)( Signal* aSignal );
            double* GetFIRFilter(int nskips);
            int GetNFilterBins(double* filterarray);
            double GetFIRSample(double* filterarray, int nfilterbins, double dtfilter, unsigned channel, unsigned patch, double AcquisitionRate);

            void InitializeBuffers(unsigned filterbuffersize, unsigned fieldbuffersize);
            void FillBuffers(Signal* aSignal, double FieldAmplitude, double FieldPhase, double LOPhase, unsigned index, unsigned channel, unsigned patch, unsigned dtauConvolutionTime);
            void PopBuffers(unsigned channel, unsigned patch);
            void CleanupBuffers();

            double* filterarray;
            double fRF_frequency;
            double fLO_frequency;
            double fAmplitude;
            double fFilter_resolution;
            std::string gfilter_filename;
            unsigned fFieldBufferSize;
            unsigned fFieldBufferMargin;
            unsigned fNPatches; // placeholder for buffer dimensioning.


            std::vector<std::deque<double>> EFieldBuffer;
            std::vector<std::deque<double>> EPhaseBuffer;
            std::vector<std::deque<double>> EAmplitudeBuffer;
            std::vector<std::deque<double>> EFrequencyBuffer;
            std::vector<std::deque<double>> LOPhaseBuffer;
            std::vector<std::deque<unsigned>> IndexBuffer;
            std::vector<std::deque<double>> PatchFIRBuffer;
            std::vector<std::deque<unsigned>> ConvolutionTimeBuffer;
    };

} /* namespace locust */

#endif /* LMCTestFIRFilterGENERATOR_HH_ */

