/*
 * LMCUtility.hh
 *
 *  Created on: Jan. 24, 2020
 *      Author: pslocum
 */

#ifndef LMCUTILITY_HH_
#define LMCUTILITY_HH_

#include "param.hh"

#include "logger.hh"

#include <vector>

namespace locust
{
 /*!
 @class Utility
 @author P. Slocum
 @brief Base class to characterize Utility selection
 @details
 Available configuration options:
 No input parameters
 */


    class Utility
    {

        public:
            Utility();
            virtual ~Utility();
    };


} /* namespace locust */

#endif /* LMCUTILITY_HH_ */
