/*
 * LMCPowerCombinerParent.hh
 *
 *  Created on: Feb 25, 2020
 *      Author: pslocum
 */

#ifndef LMCPOWERCOMBINERPARENT_HH_
#define LMCPOWERCOMBINERPARENT_HH_

namespace locust
{
 /*!
 @class LMCPowerCombinerParent
 @author P. Slocum
 @brief Base class to characterize power combiners
 @details
 Available configuration options:
 No input parameters
 */
    class PowerCombinerParent
    {

        public:
            PowerCombinerParent();
            virtual ~PowerCombinerParent();

};


} /* namespace locust */

#endif
