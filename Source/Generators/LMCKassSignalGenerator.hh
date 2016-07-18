/*
 * LMCKassSignalGenerator.hh
 *
 *  Created on: Mar 12, 2014
 *      Author: nsoblath
 */

#ifndef LMCKASSSIGNALGENERATOR_HH_
#define LMCKASSSIGNALGENERATOR_HH_

#define PI 3.1415926
//#define LO_FREQUENCY 26.8730e9 // Hz  18 keV electrons in harmonic trap.
#define LO_FREQUENCY 26.2757e9 // Hz  30 keV electrons in harmonic trap.
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
            bool DoGenerate( Signal* aSignal ) const;
            void* DriveAntenna(unsigned index, Signal* aSignal, double* ImaginarySignal) const;
            void* FilterNegativeFrequencies(Signal* aSignal, double* ImaginarySignal) const;
            double ModeExcitation() const;
            double AverageModeExcitation() const;
            double FakeModeExcitation() const;
            double* EyWR42Array() const;
            double* ScaleArray(double *array, double factor) const;
            double IntEyWR42ArraySqdA(double *EyArray1, double dimx, double dimy) const;



    };

} /* namespace locust */

#endif /* LMCKASSSIGNALGENERATOR_HH_ */
