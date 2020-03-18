/*
 * LMCUniformDistribution.cc
 *
 *  Created on: Mar 12, 2020
 *      Author: nbuzinsky
 */

#include "LMCGaussianDistribution.hh"
#include <iostream>

namespace locust
{

    GaussianDistribution::GaussianDistribution(const scarab::param_node &aParam) :
        mean( 0. ),
        std_dev( 1. )
    {
        if(aParam.has("mean"))
            mean = aParam.get_value< double >( "mean", mean );

        if(aParam.has("std-dev"))
            std_dev = aParam.get_value< double >( "std-dev", std_dev );

        distribution = std::normal_distribution<double>(mean, std_dev);
    }

    double GaussianDistribution::Generate()
    {
        return distribution(fRNEngine);
    }

} /* namespace locust */
