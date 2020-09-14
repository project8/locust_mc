/*
 * LMCTestSignalGenerator.hh
 *
 *  Created on: Jan 14 2015
 *      Author: plslocum and nsoblath
 */

#ifndef LMCTESTSIGNALGENERATOR_HH_
#define LMCTESTSIGNALGENERATOR_HH_

#include "LMCGenerator.hh"
#include "LMCRunLengthCalculator.hh"

#ifdef ROOT_FOUND
    #include "LMCRootGraphWriter.hh"
    #include "LMCRootHistoWriter.hh"
#endif

namespace scarab
{
  class param_node;
}

namespace locust
{
  class Digitizer;

    /*!
     @class TestSignalGenerator
     @author P. L. Slocum

     @brief Add Sine Wave to the signal.

     @details
     Can operate in time or frequency space

     Configuration name: "test-signal"

     Available configuration options:
     - "rf-frequency": double -- Frequency of the sine wave.
     - "lo-frequency": double -- Frequency of the local oscillator.
     - "am-frequency": double -- Frequency of AM modulation, if any.
     - "am-mod-select":  int -- Select:
      	  no AM modulation (default, 0),
      	  AM modulation (1), or
      	  AM modulation followed by demodulation (2)
      	  Following discussion at http://web.mit.edu/6.02/www/s2012/handouts/14.pdf,
      	  suggest rf-frequency = 41.0e6, lo-frequency = 40.0e6, am-frequency = 50.e6., am-mod-select = 1 or 2
      	  This allows f_(rf-lo) << f_(am) for cleaner demodulation due to band limiting.
      	  Look for loss of 2 in total power on AM modulation or demodulation as in reference above.
     - "amplitude": double -- Amplitude of the sine wave.
     - "mixing-product": bool -- If true, use product of RF and LO for downmixing.  If false,
     	 calculate lower sideband at RF-LO and replace mixing product with it.  The latter is best used
     	 for sinusoidal RF signals.  Toggling this parameter is useful to see whether the upper sideband
     	 RF+LO is downmixing into the measurement band, or not.
     - "write-root-histo": bool -- Flag to determine whether an example Root histogram is generated
         and written to file.
     - "write-root-graph": bool -- Flag to determine whether an example Root TGraph is generated
         and written to file.
     - "root-filename": string -- name of output Root file.
     - "domain": string -- Determines whether the sinusoidal test signal is generated in the time 
     	 or frequency domain
    
     Available options: "time" and "freq" [default]

    */
    class TestSignalGenerator : public Generator
    {
        public:
            TestSignalGenerator( const std::string& aName = "test-signal" );
            virtual ~TestSignalGenerator();

            bool Configure( const scarab::param_node& aNode );


            void Accept( GeneratorVisitor* aVisitor ) const;

            double GetRFFrequency() const;
            void SetRFFrequency( double aFrequency );

            double GetLOFrequency() const;
            void SetLOFrequency( double aFrequency );

            double GetAMModFrequency() const;
            void SetAMModFrequency( double aFrequency );

            int GetAMModSelect() const;
            void SetAMModSelect( int aModSelection );

            bool GetMixingProduct() const;
            void SetMixingProduct( bool aMixingProduct );

            double GetAmplitude() const;
            void SetAmplitude( double aAmplitude );

            Signal::State GetDomain() const;
            void SetDomain( Signal::State aDomain );

            bool WriteRootHisto();
            bool WriteRootGraph();



        private:
            bool DoGenerate( Signal* aSignal );

            bool DoGenerateTime( Signal* aSignal );
            bool DoGenerateFreq( Signal* aSignal );

            bool (TestSignalGenerator::*fDoGenerateFunc)( Signal* aSignal );

            double fRF_frequency;
            double fLO_frequency;
            double fAMModfrequency;
            int fAMModSelect;
            double fAmplitude;
            bool fMixingProduct;
            bool fWriteRootHisto;
            bool fWriteRootGraph;
            std::string fRootFilename;

            
    };

} /* namespace locust */

#endif /* LMCTestSignalGENERATOR_HH_ */

