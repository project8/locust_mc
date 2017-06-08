/*
 * LMCKassSignalGenerator.hh
 *
 *  Created on: Mar 12, 2014
 *      Author: nsoblath
 */

#ifndef LMCKASSSIGNALGENERATOR_HH_
#define LMCKASSSIGNALGENERATOR_HH_

#define LO_FREQUENCY 0.
//#define LO_FREQUENCY 26.8730e9 // Hz  18 keV electrons in harmonic trap.
//#define LO_FREQUENCY 26.2757e9 // Hz  30 keV electrons in harmonic trap, pitch 87-90 spans 85-50 MHz in baseband. 
//#define LO_FREQUENCY 26.3057e9 // Hz  30 keV electrons in harmonic trap, pitch 86-90 spans 81-20 MHz in baseband
////#define LO_FREQUENCY 26.2524e9 // Hz  30.48 keV electrons in harmonic trap, pitch 90, sits at 50 MHz in baseband
//#define LO_FREQUENCY 24.8100e9 // Hz, 30.48 keV electron in harmonic trap, pitch 90, main field 0.9454, sits at 50 MHz in baseband.
//#define LO_FREQUENCY 20.9688e9 // Hz, 30.48 keV electron in harmonic trap, pitch 90, main field 0.8 T, sits at 50 MHz in baseband.

//#define LO_FREQUENCY 25.1159e9 // Hz, 30.48 keV electron in harmonic trap, pitch 90, main field 0.9583 as per Luiz.  Sits at 84.9 MHz in baseband.

//#define LO_FREQUENCY 25.1010e9 // Hz, 30.48 keV electron in harmonic trap, pitch 90, main field 0.9583 as per Luiz.  Sits at 99.8 MHz in baseband.               




//#define LO_FREQUENCY 25.7537e9 // Hz, 17.83 keV electron in harmonic trap, pitch 90, main field 0.9583 as per Luiz.  Sits at 50 MHz in baseband.               



//#define LO_FREQUENCY 25.14e9 // Hz, same frequency as RSA
//#define LO_FREQUENCY 25.7737e9 // Hz, mix central 17.83 keV down to 30MHz.
//#define LO_FREQUENCY 26.3550e9 // Hz, find upper sideband of 17.83 keV

//#define LO_FREQUENCY 22.2897e9 // Hz, 30.48 keV electron in harmonic trap, pitch 90, main field 0.85 T.  Sits at 50 MHz in baseband.
//#define LO_FREQUENCY 27.0062e9  // Hz 18 keV electrons in bathtub trap.
//#define LO_FREQUENCY 26.4061e9  // Hz 30 keV electrons in bathtub trap.
//#define LO_FREQUENCY 26.3757e9 // Hz  30 keV electrons in harmonic trap, mixed to -50 MHz.


#include "LMCGenerator.hh"

namespace locust
{

    /*!
     @class KassSignalGenerator
     @author N. S. Oblath

     @brief

     @details
     Operates in time space

     Configuration name: "kass-signal"

     Available configuration options:
     - "param-name": type -- Description

    */
    class KassSignalGenerator : public Generator
    {
        public:

            KassSignalGenerator( const std::string& aName = "kass-signal" );
            virtual ~KassSignalGenerator();

            bool Configure( const ParamNode* aNode );

            void Accept( GeneratorVisitor* aVisitor ) const;

        private:
      double fLO_Frequency;  // typically defined by a parameter in json file.
            bool DoGenerate( Signal* aSignal ) const;
            void* DriveAntenna(int PreEventCounter, unsigned index, Signal* aSignal, double* ImaginarySignal) const;
            void* FilterNegativeFrequencies(Signal* aSignal, double* ImaginarySignal) const;
            double TE11ModeExcitation() const;



    };

} /* namespace locust */

#endif /* LMCKASSSIGNALGENERATOR_HH_ */
