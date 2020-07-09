/*
 * LMCBaseDistribution.cc
 *
 *  Created on: Mar 10, 2020
 *      Author: nbuzinsky
 */

#include "LMCBaseDistribution.hh"

namespace locust
{
    BaseDistribution::BaseDistribution()
    {
    }

    BaseDistribution::~BaseDistribution()
    {
    }

    void BaseDistribution::SetRandomEngine(const std::shared_ptr<std::default_random_engine> &aEngine)
    {
        fRNEngine = aEngine;
    }

} /* namespace locust */
