/*
 * LMCPowerCombiner.cc
 *
 *  Created on: Jan 27, 2019
 *      Author: penny
 */

#include "LMCPowerCombiner.hh"


namespace locust
{

    PowerCombiner::PowerCombiner()
    {
    }


    PowerCombiner::~PowerCombiner()
    {
    }



    double PowerCombiner::GetVoltageDamping(int njunctions)
    {
    	// test:  clobber voltage with 0.6 at each junction.
    return njunctions*0.6;
    }

    double PowerCombiner::GetLinePhaseCorr(unsigned z_index, double DopplerFrequency)
    {
        double D = 0.007711*0.95; // m.  18.0 keV 90 degree electron, lambda in kapton.
        double c_n = LMCConst::C()/1.5;  // speed of light in Kapton.
        double lambda = c_n/(DopplerFrequency/(2.*LMCConst::Pi()));
        double dphi = 2.*LMCConst::Pi() * D * z_index / lambda;
        return dphi;
    }




}


