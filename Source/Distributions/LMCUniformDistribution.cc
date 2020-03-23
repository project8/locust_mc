/*
 * LMCUniformDistribution.cc
 *
 *  Created on: Mar 10, 2020
 *      Author: nbuzinsky
 */

#include "LMCUniformDistribution.hh"
#include <iostream>

namespace locust
{

    UniformDistribution::UniformDistribution(const scarab::param_node &aParam) :
        fMinValue( 0. ),
        fMaxValue( 1. )
    {
        if(aParam.has("min-value"))
            fMinValue = aParam.get_value< double >( "min-value", fMinValue );

        if(aParam.has("max-value"))
            fMaxValue = aParam.get_value< double >( "max-value", fMaxValue );

        fDistribution = std::uniform_real_distribution<double>(fMinValue, fMaxValue);
    }

    double UniformDistribution::Generate()
    {
        return fDistribution(fRNEngine);
    }

} /* namespace locust */
