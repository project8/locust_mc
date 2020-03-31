/*
 * LMCUniformDistribution.cc
 *
 *  Created on: Mar 12, 2020
 *      Author: nbuzinsky
 */

#include "LMCGaussianDistribution.hh"

namespace locust
{
    LOGGER( lmclog, "LMCGaussianDistribution" );

    GaussianDistribution::GaussianDistribution(const scarab::param_node &aParam) :
        fMean( 0. ),
        fStdDev( 1. )
    {
        if(aParam.has("mean"))
            fMean = aParam.get_value< double >( "mean", fMean );

        if(aParam.has("std-dev"))
            fStdDev = aParam.get_value< double >( "std-dev", fStdDev );

        fDistribution = std::normal_distribution<double>(fMean, fStdDev);
        LDEBUG( lmclog, "Created gaussian distribution. mean: " <<fMean<<" std-dev: "<<fStdDev);
    }

    GaussianDistribution::GaussianDistribution(const double &aMean, const double &aStdDev) :
        fMean(aMean),
        fStdDev(aStdDev),
        fDistribution(aMean, aStdDev)
    {
    }

    double GaussianDistribution::Generate()
    {
        return fDistribution(*fRNEngine);
    }

} /* namespace locust */
