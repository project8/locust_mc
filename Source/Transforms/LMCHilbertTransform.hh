/*
 * LMCHilbertTransform.hh
 *
 *  Created on: May 20, 2019
 *      Author: pslocum
 */

#ifndef LMCHILBERTTRANSFORM_HH_
#define LMCHILBERTTRANSFORM_HH_

#include <fftw3.h>
#include <math.h>
#include <deque>
#include "LMCConst.hh"
#include "param.hh"


namespace locust
{
 /*!
 @class HilbertTransform
 @author P. Slocum
 @brief Class to handle Hilbert transforms for deriving mag, phase, and mean of arbitrary incident real fields on arrival at receiver.
 @details
 Available configuration options:  none yet.
 No input parameters
 */
    class HilbertTransform
    {

        public:
            HilbertTransform();
            virtual ~HilbertTransform();
            bool Configure( const scarab::param_node& aNode);
            double* GetMagPhaseMean(std::deque<double> FieldBuffer, std::deque<double> FrequencyBuffer);
            void SetBufferSize( int aBufferSize );
            void SetBufferMargin( int aBufferMargin );
            int GetBufferSize();
            int GetBufferMargin();

        private:

            fftw_complex* Transform(std::deque<double> FieldBuffer);
            double* GetFrequencyData(std::deque<double> FrequencyBuffer);
            double GetMean( std::deque<double> FieldBuffer );
            double GetMean( fftw_complex* array, int IQ, int size );
            double* GetSpan( fftw_complex* array, int IQ, int size );
            double GetPhase( double VI, double VQ, double VMean);
            double QuadrantCorrection( double VI, double HilbertPhase, double HilbertMean );
            int fbufferMargin;
            int fbufferSize;

    };

} /* namespace locust */

#endif /* LMCHILBERTTRANSFORM_HH_ */
