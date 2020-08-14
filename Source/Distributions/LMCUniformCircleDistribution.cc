/*
 * LMCUniformCircleDistribution.cc
 *
 *  Created on: Aug 12, 2020
 *      Author: nbuzinsky
 */

#include "LMCUniformCircleDistribution.hh"

namespace locust
{
    LOGGER( lmclog, "LMCUniformCircleDistribution" );

    UniformCircleDistribution::UniformCircleDistribution(const scarab::param_node &aParam) :
        fRadius( 0. )
    {
        if(aParam.has("radius"))
            fRadius = aParam.get_value< double >( "radius", fRadius );

        fDistribution = std::uniform_real_distribution<double>(0, 1);
        LDEBUG( lmclog, "Created uniform distribution. radius: " <<fRadius);
    }

    UniformCircleDistribution::UniformCircleDistribution(const double &aRadius) :
        fRadius(aRadius),
        fDistribution(0, 1)
    {
    }

    double UniformCircleDistribution::Generate()
    {
        return fRadius * sqrt( fDistribution(*fRNEngine));
    }

} /* namespace locust */
