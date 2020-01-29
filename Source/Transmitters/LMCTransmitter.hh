/*
 * LMCTransmitter.hh
 *
 *  Created on: Jan. 24, 2020
 *      Author: pslocum
 */

#ifndef LMCTRANSMITTER_HH_
#define LMCTRANSMITTER_HH_

#include "LMCThreeVector.hh"
#include "LMCPatchAntenna.hh"
#include "LMCSlotAntenna.hh"


namespace locust
{
 /*!
 @class Transmitter
 @author P. Slocum
 @brief Base class to characterize Transmitter selection
 @details
 Available configuration options:
 No input parameters
 */


    class Transmitter
    {

        public:
            Transmitter();
            virtual ~Transmitter();
            virtual void TxSayHello();

            virtual double* GetEFieldCoPol(Receiver* currentElement, int channelIndex, int zIndex, double elementSpacing, int nElementsPerStrip, double dt) {};

            virtual bool IsKassiopeia() {return false;};

        private:


};


} /* namespace locust */

#endif /* LMCCHANNEL_HH_ */
