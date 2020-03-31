/*
 * LMCDiracDistribution.cc
 *
 *  Created on: Mar 24, 2020
 *      Author: nbuzinsky
 */

#include "LMCDiracDistribution.hh"

namespace locust
{
    LOGGER( lmclog, "LMCDiracDistribution" );

    DiracDistribution::DiracDistribution(const scarab::param_node &aParam) :
        fValue( 0. )
    {
        if(aParam.has("value"))
            fValue = aParam.get_value< double >( "value", fValue );

        LDEBUG( lmclog, "Created dirac distribution. value: " << fValue);
    }

    DiracDistribution::DiracDistribution(const double &aValue) :
        fValue(aValue)
    {
    }

    double DiracDistribution::Generate()
    {
        return fValue;
    }

} /* namespace locust */
