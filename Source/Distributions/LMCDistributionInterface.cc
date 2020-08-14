/*
 * LMCDistributionInterface.cc
 *
 *  Created on: Mar 11, 2020
 *      Author: nbuzinsky
 */

#include "LMCDistributionInterface.hh"

namespace locust
{
    LOGGER( lmclog, "LMCDistributionInterface" );

    DistributionInterface::DistributionInterface() :
        fSeed(0)
    {
        fRNEngine = std::make_shared<std::default_random_engine>();
        fRNEngine->seed(fSeed);
    }

    DistributionInterface::~DistributionInterface()
    {
    }

    void DistributionInterface::SetSeed(const unsigned &aSeed)
    {
        fSeed = aSeed;
        fRNEngine->seed(fSeed);
    }

    std::shared_ptr< BaseDistribution> DistributionInterface::get_dist(const scarab::param_node &aParam)
    {
        std::string dist_name;
        if(aParam.has("name"))
        {
            dist_name = aParam.get_value< std::string >( "name", dist_name );
        }
        else
        {
            LERROR(lmclog,"Distribution name undefined!");
        }

        if(dist_name == "dirac" || dist_name == "fixed" )
            fDistributionList.push_back( std::make_shared< DiracDistribution >(aParam) );

        else if(dist_name == "exponential")
            fDistributionList.push_back( std::make_shared< ExponentialDistribution >(aParam) );

        else if(dist_name == "gaussian")
            fDistributionList.push_back( std::make_shared< GaussianDistribution >(aParam) );

        else if(dist_name == "lorentzian" || dist_name == "cauchy")
            fDistributionList.push_back( std::make_shared< LorentzianDistribution >(aParam) );

        else if(dist_name == "rudd")
            fDistributionList.push_back( std::make_shared< RuddDistribution >(aParam) );

        else if(dist_name == "uniform")
            fDistributionList.push_back( std::make_shared< UniformDistribution >(aParam) );

        else if(dist_name == "uniform-circle")
            fDistributionList.push_back( std::make_shared< UniformCircleDistribution >(aParam) );

        else if(dist_name == "kr-complex-line")
            fDistributionList.push_back( std::make_shared< KrComplexLineDistribution >(aParam) );

        else
        {
            LERROR(lmclog,"Distribution " << dist_name << " not found in list of possibilities! Enter a valid distribution name!");
        }

        fDistributionList.back()->SetRandomEngine(fRNEngine);

        return fDistributionList.back();  
    }

} /* namespace locust */
