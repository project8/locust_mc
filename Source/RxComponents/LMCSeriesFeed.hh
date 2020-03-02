/*
 * LMCSeriesFeed.hh
 *
 *  Created on: Feb 28, 2020
 *      Author: pslocum
 */

#ifndef LMCSERIESFEED_HH_
#define LMCSERIESFEED_HH_

#include "LMCPowerCombinerParent.hh"
#include "param.hh"
#include "logger.hh"
#include <iostream>

namespace locust
{
 /*!
 @class SeriesFeed
 @author P. Slocum
 @brief Derived class describing a series feed.
 @details
 Available configuration options:
 No input parameters
 */


    class SeriesFeed: public PowerCombinerParent
    {

        public:
            SeriesFeed();
            virtual ~SeriesFeed();
            virtual bool Configure( const scarab::param_node& aNode );
        	virtual bool SetVoltageDampingFactors();



        private:





    };


} /* namespace locust */

#endif /* LMCSERIESFEED_HH_ */
