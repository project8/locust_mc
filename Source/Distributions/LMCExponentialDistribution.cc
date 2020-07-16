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
        fTau( 1. ),
        fShift( 0. )
    {
        if(aParam.has("tau"))
            fTau = aParam.get_value< double >( "tau", fTau );
        if(aParam.has("shift"))
            fShift = aParam.get_value< double >( "shift", fShift );

        fDistribution = std::exponential_distribution<double>(1. / fTau);
        LDEBUG( lmclog, "Created exponential distribution. tau: " <<fTau<<" shift: "<<fShift);
    }

    ExponentialDistribution::ExponentialDistribution(const double &aTau, const double &aShift) :
        fTau( aTau ),
        fShift( aShift ),
        fDistribution( aTau )
    {
    }

    double ExponentialDistribution::Generate()
    {
        return fDistribution(*fRNEngine) + fShift;
    }

} /* namespace locust */
