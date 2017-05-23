/*
 * LMCKassSignalGenerator.hh
 *
 *  Created on: Mar 12, 2014
 *      Author: nsoblath
 */

#ifndef LMCKASSSIGNALGENERATOR_HH_
#define LMCKASSSIGNALGENERATOR_HH_

#define LO_FREQUENCY 0.  // This will be overwritten by parameter from json file.

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
            double ModeExcitation() const;
            double AverageModeExcitation() const;
            double FakeModeExcitation() const;
            double* EyWR42Array() const;
            double* ScaleArray(double *array, double factor) const;
            double IntEyWR42ArraySqdA(double *EyArray1, double dimx, double dimy) const;
            double Phase2TE11ModeExcitation() const;



    };

} /* namespace locust */

#endif /* LMCKASSSIGNALGENERATOR_HH_ */
