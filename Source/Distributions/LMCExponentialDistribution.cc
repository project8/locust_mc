/*
 * LMCExponentialDistribution.cc
 *
 *  Created on: Mar 23, 2020
 *      Author: nbuzinsky
 */

#include "LMCExponentialDistribution.hh"

namespace locust
{
    LOGGER( lmclog, "LMCExponentialDistribution" );

    ExponentialDistribution::ExponentialDistribution(const scarab::param_node &aParam) :
        fLambda( 1. ),
        fShift( 0. )
    {
        if(aParam.has("lambda"))
            fLambda = aParam.get_value< double >( "lambda", fLambda );
        if(aParam.has("shift"))
            fShift = aParam.get_value< double >( "shift", fShift );

        fDistribution = std::exponential_distribution<double>(fLambda);
        LDEBUG( lmclog, "Created exponential distribution. lambda: " <<fLambda<<" shift: "<<fShift);
    }

    ExponentialDistribution::ExponentialDistribution(const double &aLambda, const double &aShift) :
        fLambda( aLambda ),
        fShift( aShift ),
        fDistribution( aLambda )
    {
    }

    double ExponentialDistribution::Generate()
    {
        return fDistribution(*fRNEngine) + fShift;
    }

} /* namespace locust */
