/*
 * LMCDipoleAntenna.hh
 *
 *  Created on: Dec 19, 2019
 *      Author: pslocum
 */

#ifndef LMCDIPOLEANTENNA_HH_
#define LMCDIPOLEANTENNA_HH_

#include "LMCTransmitterHardware.hh"


namespace locust
{
 /*!
 @class DipoleAntenna
 @author P. Slocum
 @brief Derived class describing a dipole antenna transmitter.
 @details
 Available configuration options:
 No input parameters
 */
    class DipoleAntenna: public TransmitterHardware
    {

        public:
            DipoleAntenna();
            virtual ~DipoleAntenna();

            virtual void TxHardwareSayHello();

            virtual double GetPatternFactor(LMCThreeVector pointOfInterest);

        private:

            LMCThreeVector fMomentVector;


    };


} /* namespace locust */

#endif /* LMCDIPOLEANTENNA_HH_ */
