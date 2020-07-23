/*
 * LMCAntennaElementPositioner.hh
 *
 *  Created on: Feb 25, 2020
 *      Author: pslocum
 */

#ifndef LMCANTENNAELEMENTPOSITIONER_HH_
#define LMCANTENNAELEMENTPOSITIONER_HH_
#include "param.hh"
#include "LMCException.hh"


namespace locust
{
 /*!
 @class LMCAntennaElementPositioner
 @author P. Slocum
 @brief Base class to characterize power combiners
 @details
 Available configuration options:
 No input parameters
 */
    class AntennaElementPositioner
    {

        public:
            AntennaElementPositioner();
            virtual ~AntennaElementPositioner();
            virtual bool Configure( const scarab::param_node& aNode );

};


} /* namespace locust */

#endif
