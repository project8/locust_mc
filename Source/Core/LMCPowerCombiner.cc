/*
 * LMCPowerCombiner.cc
 *
 *  Created on: Jan 27, 2019
 *      Author: penny
 */

#include "LMCPowerCombiner.hh"
#include <iostream>


namespace locust
{

    PowerCombiner::PowerCombiner()
    {
    }


    PowerCombiner::~PowerCombiner()
    {
    }

    double PowerCombiner::GetCorporateVoltageDamping()
    {
    	// currently hard-coded for one-quarter power combining

    	double dampingfactor = 0.;
    	dampingfactor = 0.38*0.66; // the patch and the amplifier only.
    	return dampingfactor;
    }

    // **NOTE FOR THE FOLLOWING ROUTINES**
    // z_index is the index of the current patch, starting at zero.
    // njunctions is the number of junctions the signal has to pass through between this patch and the amp.
    // for ex, if I am on the third patch to the right of the amp, njunctions = 2.

    double PowerCombiner::GetSeriesPhaseDelay(unsigned z_index, double DopplerFrequency, double PatchSpacing)
    {
        int njunctions = z_index;
        double D = PatchSpacing; // m.  18.6 keV 90 degree electron, lambda in kapton.
        double c_n = LMCConst::C()/1.704;  // speed of light in Kapton.
        double lambda = c_n/DopplerFrequency;
        double dphi = 2.*LMCConst::Pi() * D * njunctions / lambda;
        return dphi;
    }

    double PowerCombiner::GetSeriesVoltageDamping(unsigned z_index)
    {
    	// currently hard-coded one quarter power combining
    	int njunctions = z_index;
    	double dampingfactor = 0.;
    	dampingfactor = 0.38*0.66*pow(0.87, njunctions);
		return dampingfactor;
    }

    double PowerCombiner::GetCenterFedPhaseDelay(int NPatchesPerStrip, unsigned z_index, double DopplerFrequency, double PatchSpacing)
      {
    	int njunctions = fabs((double)z_index - (double)NPatchesPerStrip/2.) - 1;
    	if (z_index >= NPatchesPerStrip/2) njunctions += 1; // compensate for patches to the right of amp.
        // printf("z_index is %d and njunctions is %d\n", z_index, njunctions); getchar();
        double D = PatchSpacing; // m.  18.0 keV 90 degree electron, lambda in kapton.
        double c_n = LMCConst::C()/1.704;  // speed of light in Kapton.
        double lambda = c_n/DopplerFrequency;
        double dphi = 2.*LMCConst::Pi() * D * njunctions / lambda;
        return dphi;
      }


    double PowerCombiner::GetOneQuarterVoltageDamping(int NPatchesPerStrip, unsigned z_index)
      {
          int njunctions = fabs((double)z_index - (double)NPatchesPerStrip/2.) - 1;
          if (z_index >= NPatchesPerStrip/2) njunctions += 1; // compensate for patches to the right of amp.
          double dampingfactor = 0.;
	      dampingfactor = 0.38*pow(0.87, njunctions)*0.66; // patch loss * junction loss * amplifier loss
         // printf("dampingfactor is %g\n", dampingfactor); getchar();
	      return dampingfactor;
      }

    double PowerCombiner::GetSevenEighthsVoltageDamping(int NPatchesPerStrip, unsigned z_index)
          {
    		int njunctions = fabs((double)z_index - (double)NPatchesPerStrip/2.) - 1;
    	    if (z_index >= NPatchesPerStrip/2) njunctions += 1; // compensate for patches to the right of amp.
    	    double dampingfactor = 0.;
    		dampingfactor = 0.6*pow(0.75, njunctions)*0.66; // patch loss * junction loss * amplifier loss
         //   printf("dampingfactor is %g\n", dampingfactor); getchar();
    	    return dampingfactor;
          }

    double PowerCombiner::GetNineSixteenthsVoltageDamping(int NPatchesPerStrip, unsigned z_index)
    {
    	int njunctions = fabs((double)z_index - (double)NPatchesPerStrip/2.) - 1;
    	if (z_index >= NPatchesPerStrip/2) njunctions += 1; // compensate for patches to the right of amp.
        double dampingfactor = 0.;
    	dampingfactor = 0.52*pow(0.8, njunctions)*0.66; // patch loss * junction loss * amplifier loss
    	//   printf("dampingfactor is %g\n", dampingfactor); getchar();
    	return dampingfactor;
    }

}


