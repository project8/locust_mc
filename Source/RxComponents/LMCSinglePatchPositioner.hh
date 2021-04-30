/*
 * LMCSinglePatchPositioner.hh
 *
 *  Created on: Feb 28, 2020
 *      Author: pslocum
 */

#ifndef LMCSINGLEPATCHPOSITIONER_HH_
#define LMCSINGLEPATCHPOSITIONER_HH_

#include "LMCAntennaElementPositioner.hh"
#include "param.hh"
#include "logger.hh"
#include <iostream>

namespace locust
{
 /*!
 @class SinglePatchPositioner
 @author P. Slocum
 @brief Derived class describing a single patch.
 @details
 Available configuration options:
 No input parameters
 */


    class SinglePatchPositioner: public AntennaElementPositioner
    {

        public:
            SinglePatchPositioner();
            virtual ~SinglePatchPositioner();
            virtual bool Configure( const scarab::param_node& aNode );

            virtual double GetPositionZ(double zShiftArray, int channelIndex, int nChannels,
             		int nSubarrays, int nReceivers, double elementSpacingZ, int receiverIndex);

    };


} /* namespace locust */

#endif /* LMCSINGLEPATCHPOSITIONER_HH_ */
