/*
 * LMCField.hh
 *
 *  Created on: Jun. 4, 2021
 *      Author: pslocum
 */

#ifndef LMCFIELD_HH_
#define LMCFIELD_HH_

#include "param.hh"

#include "logger.hh"

#include <vector>

namespace locust
{
 /*!
 @class Field
 @author P. Slocum
 @brief Base class to characterize Field selection
 @details
 Available configuration options:
 No input parameters
 */


    class Field
    {

        public:
            Field();
            virtual ~Field();
            virtual void SayHello();

            virtual bool Configure( const scarab::param_node& ){};
    };


}; /* namespace locust */

#endif /* LMCFIELD_HH_ */
