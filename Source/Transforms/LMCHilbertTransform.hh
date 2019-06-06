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
            double* GetMagPhaseMean(std::deque<double> FieldBuffer, std::deque<double> FrequencyBuffer, unsigned edge_margin, double AcquisitionRate);

        private:

            fftw_complex* Transform(std::deque<double> FieldBuffer);
            double* GetFrequencyData(std::deque<double> FrequencyBuffer);
            double GetMean( std::deque<double> FieldBuffer );
            double QuadrantCorrection( std::deque<double> FieldBuffer, double HilbertPhase, double HilbertMean );
    };

} /* namespace locust */

#endif /* LMCHILBERTTRANSFORM_HH_ */
