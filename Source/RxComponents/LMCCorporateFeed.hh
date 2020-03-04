/*
 * LMCCorporateFeed.hh
 *
 *  Created on: Feb 28, 2020
 *      Author: pslocum
 */

#ifndef LMCCORPORATEFEED_HH_
#define LMCCORPORATEFEED_HH_

#include "LMCPowerCombiner.hh"
#include "param.hh"
#include "logger.hh"
#include <iostream>

namespace locust
{
 /*!
 @class CorporateFeed
 @author P. Slocum
 @brief Derived class describing a single patch.
 @details
 Available configuration options:
 No input parameters
 */


    class CorporateFeed: public PowerCombiner
    {

        public:
            CorporateFeed();
            virtual ~CorporateFeed();
            virtual bool Configure( const scarab::param_node& aNode );

        	virtual bool SetVoltageDampingFactors();



    };


} /* namespace locust */

#endif /* LMCCORPORATEFEED_HH_ */
