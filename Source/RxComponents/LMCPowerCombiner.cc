/*
 * LMCPowerCombiner.cc
 *
 *  Created on: Jan 27, 2019
 *      Author: p. l. slocum
 */

#include "LMCPowerCombiner.hh"
#include <iostream>


namespace locust
{

  PowerCombiner::PowerCombiner():
  nPatchesPerStrip( 0 ),
  junctionLoss( 0 ),
  patchLoss( 0 ),
  amplifierLoss( 0 ),
  endPatchLoss( 0 ),
  dampingFactors( 0 )
  {
  }


  PowerCombiner::~PowerCombiner()
  {
  }


  bool PowerCombiner::AddOneVoltageToStripSum(Signal* aSignal, double VoltageFIRSample, double phi_LO, unsigned z_index, unsigned sampleIndex)
  {
      VoltageFIRSample *= dampingFactors[z_index];

	  aSignal->LongSignalTimeComplex()[sampleIndex][0] += 2.*VoltageFIRSample * sin(phi_LO);
      aSignal->LongSignalTimeComplex()[sampleIndex][1] += 2.*VoltageFIRSample * cos(phi_LO);

      return true;

  }





  double PowerCombiner::GetVoltageDividerWeight(double RJunction, double R0, double RGround, int NPatchesPerStrip, unsigned z_index)
  {
	int NPAIRS = fabs((double)z_index - (double)NPatchesPerStrip/2.);
	if (z_index >= NPatchesPerStrip/2) NPAIRS += 1; // compensate for patches to the right of amp.

	std::vector<double> D = GetPartialGains(RJunction, R0, RGround, NPAIRS);  // calculate new vector of gains.
	double dampingfactor = 0.6*0.66;  // patch loss * T-junction loss.
	return dampingfactor * D[NPAIRS-1];
  }

  std::vector<double> PowerCombiner::GetResistances(double RJunction, double R0, double RGround, int NPAIRS)
  {
	int NRESISTORS = NPAIRS + 1;
	std::vector<double> R;
	R.resize(NRESISTORS);

	for (unsigned i=0; i<NPAIRS; i++)
	  {
	  R[i] = R0 + i*RJunction;
	  }
	R[NRESISTORS-1] = RGround;

	return R;
  }

  std::vector<double> PowerCombiner::GetPartialGains(double RJunction, double R0, double RGround, int NPAIRS)
  {
	int NRESISTORS = NPAIRS + 1;
	std::vector<double> D;
	D.resize(NRESISTORS);
	std::vector<double> R = GetResistances(RJunction, R0, RGround, NPAIRS);
	double Dtotal = 0.;

	for (int i=0; i<NRESISTORS; i++)
	  {
	  double ParallelResistance = GetParallelResistance(R, NRESISTORS, i);
	  D[i] = ParallelResistance/(ParallelResistance + R[i]);
	  Dtotal += D[i];
//	  printf("D[%d] is %f and Dtotal is %f\n", i, D[i], Dtotal);
	  }
//	getchar();

	return D;

  }
  double PowerCombiner::GetParallelResistance(std::vector<double> R, int NRESISTORS, int resistorindex)
  {

	double OneOverR = 0.;
	for (unsigned i=0; i<NRESISTORS; i++)
	  {
	  if (i != resistorindex) OneOverR += 1./R[i];
	  }
	return 1./OneOverR;

  }


/*
  double PowerCombiner::GetCorporateVoltageDamping()
  {
    // currently hard-coded for one-quarter power combining

    double dampingfactor = 0.;
    dampingfactor = 0.38*0.66; // the patch and the amplifier only.
    return dampingfactor;
  }

  */

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



  void PowerCombiner::SetVoltageDividerDampingFactors()
  {
	  for (unsigned z_index=0; z_index<nPatchesPerStrip; z_index++)
	  {
		  int NPAIRS = fabs((double)z_index - (double)nPatchesPerStrip/2.);
		  if (z_index >= nPatchesPerStrip/2) NPAIRS += 1; // compensate for patches to the right of amp.
		  std::vector<double> D = GetPartialGains(junctionResistance, 1.0, 10.e6, NPAIRS);  // calculate new vector of gains.
		  dampingFactors[z_index] = patchLoss*amplifierLoss * D[NPAIRS-1];  // patch loss * T-junction loss
	  }
  }


  void PowerCombiner::SetSeriesFedDampingFactors()
  {

      for (unsigned z_index=0; z_index<nPatchesPerStrip; z_index++)
      {
          int njunctions = z_index;
          dampingFactors[z_index] = patchLoss*amplifierLoss*pow(junctionLoss, njunctions);
      }

  }

/*
  double PowerCombiner::GetSeriesVoltageDamping(unsigned z_index)
  {
    // currently hard-coded one quarter power combining
    int njunctions = z_index;
    double dampingfactor = 0.;
    dampingfactor = 0.38*0.66*pow(0.87, njunctions);
    return dampingfactor;
  }
*/

  double PowerCombiner::GetCenterFedPhaseDelay(int NPatchesPerStrip, unsigned z_index, double DopplerFrequency, double PatchSpacing)
  {
    int njunctions = fabs((double)z_index - (double)NPatchesPerStrip/2.) - 1;
    if (z_index >= NPatchesPerStrip/2) njunctions += 1; // compensate for patches to the right of amp.
    double D = PatchSpacing; // m.  18.0 keV 90 degree electron, lambda in kapton.
    double c_n = LMCConst::C()/1.704;  // speed of light in Kapton.
    double lambda = c_n/DopplerFrequency;
    double dphi = 2.*LMCConst::Pi() * D * njunctions / lambda;
    return dphi;
  }

  void PowerCombiner::SetCenterFedDampingFactors()
  {

	  for (unsigned z_index=0; z_index<nPatchesPerStrip; z_index++)
	  {
	      int njunctions = fabs((double)z_index - (double)nPatchesPerStrip/2.) - 1;
	      if (z_index >= nPatchesPerStrip/2) njunctions += 1; // compensate for patches to the right of amp.
	      if (z_index == 0 || z_index == nPatchesPerStrip-1)
	      	  {
	    	  dampingFactors[z_index] = endPatchLoss*pow(junctionLoss, njunctions)*amplifierLoss;
	      	  }
	      else
	      	  {
	    	  dampingFactors[z_index] = patchLoss*pow(junctionLoss, njunctions)*amplifierLoss; // patch loss * junction loss * amplifier loss
	      	  }
	  }

  }

  void PowerCombiner::SetVoltageDampingFactors(int aPowerCombiner, int aPatchesPerStrip)
  {

	  SetNPatchesPerStrip(aPatchesPerStrip);
	  dampingFactors.resize(nPatchesPerStrip);

	  if ((aPowerCombiner == 0) || (aPowerCombiner == 2) || (aPowerCombiner == 3) || (aPowerCombiner == 4))
	  {
		  SetCenterFedDampingFactors();
	  }

	  else if (aPowerCombiner == 1)
	  {
		  SetSeriesFedDampingFactors();
	  }

	  else if (aPowerCombiner == 5)
	  {
          SetVoltageDividerDampingFactors();
	  }

  }


/*
  double PowerCombiner::GetOneQuarterVoltageDamping(int NPatchesPerStrip, unsigned z_index)
  {
    int njunctions = fabs((double)z_index - (double)NPatchesPerStrip/2.) - 1;
    if (z_index >= NPatchesPerStrip/2) njunctions += 1; // compensate for patches to the right of amp.
    double dampingfactor = 0.;
    if(z_index == 0 || z_index == NPatchesPerStrip-1)
      {
	dampingfactor = 0.95*pow(0.87, njunctions)*0.66;
      }
    else
      {
	dampingfactor = 0.38*pow(0.87, njunctions)*0.66; // patch loss * junction loss * amplifier loss
      }
    // printf("dampingfactor is %g\n", dampingfactor); getchar();
    return dampingfactor;
  }

  double PowerCombiner::GetSevenEighthsVoltageDamping(int NPatchesPerStrip, unsigned z_index)
  {
    int njunctions = fabs((double)z_index - (double)NPatchesPerStrip/2.) - 1;
    if (z_index >= NPatchesPerStrip/2) njunctions += 1; // compensate for patches to the right of amp.
    double dampingfactor = 0.;
    // adjust last patch to a 2-port junction.
    if(z_index == 0 || z_index == NPatchesPerStrip-1)
      {
  	dampingfactor = 0.95*pow(0.75, njunctions)*0.66;
      }
    else
      {
	dampingfactor = 0.6*pow(0.75, njunctions)*0.66; // patch loss * junction loss * amplifier loss
      }
    //   printf("dampingfactor is %g\n", dampingfactor); getchar();
    return dampingfactor;
  }

  double PowerCombiner::GetNineSixteenthsVoltageDamping(int NPatchesPerStrip, unsigned z_index)
  {
    int njunctions = fabs((double)z_index - (double)NPatchesPerStrip/2.) - 1;
    if (z_index >= NPatchesPerStrip/2) njunctions += 1; // compensate for patches to the right of amp.
    double dampingfactor = 0.;
    if(z_index == 0 || z_index == NPatchesPerStrip-1)
      {
	dampingfactor = 0.95*pow(0.8, njunctions)*0.66;
	  }
    else
      {
	dampingfactor = 0.52*pow(0.8, njunctions)*0.66; // patch loss * junction loss * amplifier loss
      }
    //   printf("dampingfactor is %g\n", dampingfactor); getchar();
    return dampingfactor;
  }
*/

  void PowerCombiner::SetSMatrixParameters(int powerCombiner, int aPatchesPerStrip)
  {
	  nPatchesPerStrip = aPatchesPerStrip;

	  if (powerCombiner == 0) // corporate
	  {
		  junctionLoss = 0.;
		  patchLoss = 0.6;
	      amplifierLoss = 0.66;
	      endPatchLoss = 0.6;
	  }

	  else if (powerCombiner == 1) // series
	  {
		  junctionLoss = 0.87;
		  patchLoss = 0.38;
	      amplifierLoss = 0.66;
	      endPatchLoss = 0.38;
	  }

	  else if (powerCombiner == 2) // one-quarter
	  {
		  junctionLoss = 0.87;
		  patchLoss = 0.38;
	      amplifierLoss = 0.66;
	      endPatchLoss = 0.95;
	  }

	  else if (powerCombiner == 3) // seven-eighths
	  {
		  junctionLoss = 0.75;
		  patchLoss = 0.6;
	      amplifierLoss = 0.66;
	      endPatchLoss = 0.95;
	  }

	  else if (powerCombiner == 4) // nine-sixteenths
	  {
		  junctionLoss = 0.8;
		  patchLoss = 0.52;
	      amplifierLoss = 0.66;
	      endPatchLoss = 0.95;
	  }

	  else if (powerCombiner == 5) // voltage-divider
	  {
		  junctionResistance = 0.3;
		  patchLoss = 0.6;
	      amplifierLoss = 0.66;
	  }


  }

  void PowerCombiner::SetNPatchesPerStrip(int aPatchesPerStrip)
  {
	  nPatchesPerStrip = aPatchesPerStrip;
  }
  void PowerCombiner::SetJunctionLoss(double aJunctionLoss)
  {
	  junctionLoss = aJunctionLoss;
  }
  void PowerCombiner::SetPatchLoss(double aPatchLoss)
  {
	  patchLoss = aPatchLoss;
  }
  void PowerCombiner::SetAmplifierLoss(double aAmplifierLoss)
  {
	  amplifierLoss = aAmplifierLoss;
  }
  void PowerCombiner::SetEndPatchLoss(double aEndPatchLoss)
  {
	  endPatchLoss = aEndPatchLoss;
  }


}


