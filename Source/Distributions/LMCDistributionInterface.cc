/*
 * LMCDistributionInterface.cc
 *
 *  Created on: Mar 11, 2020
 *      Author: nbuzinsky
 */

#include "LMCDistributionInterface.hh"

namespace locust
{
    DistributionInterface::DistributionInterface()
    {
    }

    DistributionInterface::~DistributionInterface()
    {
    }

    std::shared_ptr< BaseDistribution> DistributionInterface::get_dist(const scarab::param_node &aParam)
    {
        std::string dist_name;
        if(aParam.has("name"))
            dist_name = aParam.get_value< std::string >( "name", dist_name );
        //else
        //{
        //    LMCERROR(lmcdist,"Distribution " << dist_name << " not found in list of possibilities! Enter a valid distribution name!");
        //}

        if(dist_name == "gaussian")
            fDistributionList.push_back( std::make_shared< GaussianDistribution >(aParam) );

        else if(dist_name == "uniform")
            fDistributionList.push_back( std::make_shared< UniformDistribution >(aParam) );

        else if(dist_name == "kr-complex-line")
            fDistributionList.push_back( std::make_shared< KrComplexLineDistribution >(aParam) );

        //else error

        return fDistributionList.back();  
    }

} /* namespace locust */
