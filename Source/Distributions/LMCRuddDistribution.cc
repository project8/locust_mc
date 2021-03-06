/*
 * LMCRuddDistribution.cc
 *
 *  Created on: Mar 25, 2020
 *      Author: nbuzinsky
 */

#include "LMCRuddDistribution.hh"
#include "LMCConst.hh"

namespace locust
{
    LOGGER( lmclog, "LMCRuddDistribution" );

    RuddDistribution::RuddDistribution(const scarab::param_node &aParam) :
        fAlpha( 0.01 )
    {
        if(aParam.has("alpha"))
            fAlpha = aParam.get_value< double >( "alpha", fAlpha );

        fDistribution = std::uniform_real_distribution<double>(0., 1.);
        LDEBUG( lmclog, "Created rudd distribution. alpha: " <<fAlpha);
    }

    RuddDistribution::RuddDistribution(const double &aAlpha) :
        fAlpha(aAlpha),
        fDistribution(0., 1.)
    {
    }

    double RuddDistribution::Generate()
    {
        double u = fDistribution(*fRNEngine);
        return atan(sqrt(pow(fAlpha,2.) / (1. + pow(fAlpha,2.))) * tan( LMCConst::Pi() / 2. * u));
    }

} /* namespace locust */
