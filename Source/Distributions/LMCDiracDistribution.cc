/*
 * LMCDiracDistribution.cc
 *
 *  Created on: Mar 24, 2020
 *      Author: nbuzinsky
 */

#include "LMCDiracDistribution.hh"
#include <iostream>

namespace locust
{

    DiracDistribution::DiracDistribution(const scarab::param_node &aParam) :
        fValue( 0. )
    {
        if(aParam.has("value"))
            fValue = aParam.get_value< double >( "value", fValue );
    }

    double DiracDistribution::Generate()
    {
        return fValue;
    }

} /* namespace locust */
