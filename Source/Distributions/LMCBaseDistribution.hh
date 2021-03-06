/*
 * LMCBaseDistribution.hh
 *
 *  Created on: Mar 10, 2020
 *      Author: nbuzinsky
 */

#ifndef LMCBASEDISTRIBUTION_HH_
#define LMCBASEDISTRIBUTION_HH_

#include "logger.hh"
#include "param_node.hh"

#include <random>

namespace locust
{
 /*!
 @class BaseDistribution
 @author N. Buzinsky
 @brief Base class for Locust distributions, which can be called from config files, and used as parameters
 @details
 Available configuration options:
 No input parameters
 */
    //template<class T>
    class BaseDistribution
    {

        public:
            BaseDistribution();
            virtual ~BaseDistribution();

            virtual double Generate() = 0;
            void SetRandomEngine(const std::shared_ptr<std::default_random_engine> &aEngine);

            std::shared_ptr<std::default_random_engine> fRNEngine;

        private:
};


} /* namespace locust */

#endif /* LMCBASEDDISTRBUTION_HH_ */
