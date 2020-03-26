/*
 * LMCLorentzianDistribution.cc
 *
 *  Created on: Mar 24, 2020
 *      Author: nbuzinsky
 */

#include "LMCLorentzianDistribution.hh"
#include <iostream>

namespace locust
{

    LorentzianDistribution::LorentzianDistribution(const scarab::param_node &aParam) :
        fMean( 0. ),
        fFWHM( 1. )
    {
        if(aParam.has("mean"))
            fMean = aParam.get_value< double >( "mean", fMean );

        if(aParam.has("fwhm"))
            fFWHM = aParam.get_value< double >( "fwhm", fFWHM );

        fDistribution = std::cauchy_distribution<double>(fMean, fFWHM / 2.);
    }

    double LorentzianDistribution::Generate()
    {
        return fDistribution(*fRNEngine);
    }

} /* namespace locust */
