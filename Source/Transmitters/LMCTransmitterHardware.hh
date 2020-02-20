/*
 * LMCTransmitterHardware.hh
 *
 *  Created on: Feb. 18, 2020
 *      Author: pslocum
 */

#ifndef LMCTRANSMITTERHARDWARE_HH_
#define LMCTRANSMITTERHARDWARE_HH_

#include "param.hh"
#include "LMCThreeVector.hh"
#include "LMCConst.hh"
#include <vector>

namespace locust
{
 /*!
 @class TransmitterHardware
 @author P. Slocum
 @brief Base class to characterize TransmitterHardware (antenna, orientation in lab) selection
 @details
 Available configuration options:
 No input parameters
 */


    class TransmitterHardware
    {

        public:
            TransmitterHardware();
            virtual ~TransmitterHardware();
            virtual void TxHardwareSayHello();

            virtual double GetPatternFactor(LMCThreeVector pointOfInterest) {};

            int GetNAntennas();
            void SetNAntennas(int aNumber);

            virtual double GetDrivePhaseDifference() {};


        private:

            int fNAntennas;
            double fDrivePhaseDifference;



};


} /* namespace locust */

#endif
