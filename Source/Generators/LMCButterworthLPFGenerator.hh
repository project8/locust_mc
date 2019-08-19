/*
 * LMCButterworthLPFGenerator.hh
 *
 *  Created on: 29 January 2015
 *      Author: plslocum
 */

#ifndef LMCBUTTERWORTHLPFGENERATOR_HH_
#define LMCBUTTERWORTHLPFGENERATOR_HH_

#include "LMCGenerator.hh"

namespace locust
{

    /*!
     @class ButterworthLPFGenerator
     @author P. L. Slocum

     @brief Apply digital low pass filter to signal.

     @details
     This is an 8th order Butterworth low-pass filter defined by coefficients that are
     presently hard-coded in LMCButterworthLPFGenerator::SetCoefficients().  The coefficients
     are calculated in the Mathematica notebook as in
     https://github.com/project8/scripts/blob/master/slocum/RFButterworth8thOrder.nb .

     To check the output of this filter, use the following two generators:  "gaussian-noise"
     followed by "butterworth-lpf".  The attenuation of the noise power at high frequencies
     should be clearly visible.

     This generator can only be used after decimation.  It does not stop frequency aliasing.
     The filter behavior is wc = 70.e6 Hz, 100X attenuation at 95 MHz, fs=200 MHz.

     If a Butterworth filter is desired for use before decimation, then the filter coefficients
     would need to be recalculated at the faster pre-decimation sampling frequency.  The filtering
     would be done on aSignal->LongSignalTimeComplex() instead of aSignal->SignalTimeComplex().
     This would likely require a new generator, called LMCFastButterworthLPFGenerator, or similar.

     Configuration name: "butterworth-lpf"

    */

class ButterworthLPFGenerator : public Generator
    {
        public:
            ButterworthLPFGenerator( const std::string& aName = "butterworth-lpf" );
            virtual ~ButterworthLPFGenerator();

            bool SetCoefficients();
            bool Configure( const scarab::param_node& aNode );

            void Accept( GeneratorVisitor* aVisitor ) const;

        private:
            bool DoGenerate( Signal* aSignal );

            bool DoGenerateTime( Signal* aSignal );
            bool DoGenerateFreq( Signal* aSignal );
            bool (ButterworthLPFGenerator::*fDoGenerateFunc)( Signal* aSignal );

            double fN;
            double fM;
            std::vector<double> fFIRa;
            std::vector<double> fFIRb;
    };

} /* namespace locust */

#endif /* LMCBUTTERWORTHLPFGENERATOR_HH_ */
