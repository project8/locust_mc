/*
 * LMCSlotAntenna.hh
 *
 *  Created on: Dec 19, 2019
 *      Author: pslocum
 */

#ifndef LMCSLOTANTENNA_HH_
#define LMCSLOTANTENNA_HH_

#include "LMCReceiver.hh"

namespace locust
{
 /*!
 @class SlotAntenna
 @author P. Slocum
 @brief Derived class describing the slot antenna
 @details
 Available configuration options:
 No input parameters
 */
    class SlotAntenna: public Receiver
    {

        public:
            SlotAntenna();
            virtual ~SlotAntenna();


    };


} /* namespace locust */

#endif /* LMCSLOTANTENNA_HH_ */
