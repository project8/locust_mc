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
        fMean( 0. ),
        fStdDev( 1. )
    {
        if(aParam.has("mean"))
            fMean = aParam.get_value< double >( "mean", fMean );

        if(aParam.has("std-dev"))
            fStdDev = aParam.get_value< double >( "std-dev", fStdDev );

        fDistribution = std::normal_distribution<double>(fMean, fStdDev);
    }

    double GaussianDistribution::Generate()
    {
        return fDistribution(fRNEngine);
    }

} /* namespace locust */
