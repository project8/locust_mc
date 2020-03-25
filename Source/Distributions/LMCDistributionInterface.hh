/*
 * LMCDistributionInterface.hh
 *
 *  Created on: Mar 10, 2020
 *      Author: nbuzinsky
 */

#ifndef LMCDISTRIBUTIONINTERFACE_HH_
#define LMCDISTRIBUTIONINTERFACE_HH_

#include "LMCBaseDistribution.hh"

#include "LMCDiracDistribution.hh"
#include "LMCExponentialDistribution.hh"
#include "LMCGaussianDistribution.hh"
#include "LMCLorentzianDistribution.hh"
#include "LMCUniformDistribution.hh"
#include "LMCKrComplexLineDistribution.hh"

#include "param_node.hh"

#include <list>
#include <string>

namespace locust
{
 /*!
 @class DistributionInterface
 @author N. Buzinsky
 @brief Public interface for parsing/ calling LMCDistributions for random number generation
 @details
 Available configuration options:
 No input parameters
 */
    class DistributionInterface
    {

        public:
            DistributionInterface();
            virtual ~DistributionInterface();
            std::shared_ptr< BaseDistribution> get_dist(const scarab::param_node &aParam);

        private:
            std::list<std::shared_ptr< BaseDistribution >> fDistributionList;
};


} /* namespace locust */

#endif /* LMCBASEDISTRBUTION_HH_ */
