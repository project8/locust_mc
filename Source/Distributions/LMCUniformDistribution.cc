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
        min_value( 0. ),
        max_value( 1. )
    {
        if(aParam.has("min-value"))
            min_value = aParam.get_value< double >( "min-value", min_value );

        if(aParam.has("max-value"))
            max_value = aParam.get_value< double >( "max-value", min_value );

        distribution = std::uniform_real_distribution<double>(min_value, max_value);
    }

    double UniformDistribution::Generate()
    {
        return distribution(generator);
    }

} /* namespace locust */
