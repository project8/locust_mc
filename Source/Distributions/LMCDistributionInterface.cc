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
    //add second optional argument for type
    std::shared_ptr< BaseDistribution> DistributionInterface::get_dist(const std::string &dist_name)
    {
        ////put in case check
        if(dist_name == "uniform")
        {
            fDistributionList.push_back( std::make_shared< UniformDistribution >() );
        }
        //else
        //{
        //    LMCERROR(lmcdist,"Distribution " << dist_name << " not found in list of possibilities! Enter a valid distribution name!");
        //}

        return fDistributionList.back();  
    }

} /* namespace locust */
