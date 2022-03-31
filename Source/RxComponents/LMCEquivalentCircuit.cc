/*
 * LMCEquivalentCircuit.cc
 *
 *  Created on: Nov 12, 2021
 *      Author: jgaison
 */

#include "LMCEquivalentCircuit.hh"
#include <iostream>
#include <stdio.h>
#include "logger.hh"

namespace locust
{

	LOGGER( lmclog, "EquivalentCircuit" );


    EquivalentCircuit::EquivalentCircuit():
    		fEquivalentR( 1. ),
			fEquivalentL( 0.159e-6 ),
			fEquivalentC( 0.159e-12 ),
			fTFBins( 4000 ),
			fFreqRangeCenter( 1.0e9 ),
			fGeneratingTF( false )
    {}
    EquivalentCircuit::~EquivalentCircuit() {}


    bool EquivalentCircuit::Configure( const scarab::param_node& aParam )
    {
		//Update any parameters defined in the config file
		//If any of the R, L, or C parameters are updated, set fGeneratingTF to true so that TF will be generated
    	if( aParam.has( "equivalentR" ) )
    	{
    		fEquivalentR = aParam["equivalentR"]().as_double();
    		fGeneratingTF = true;
    	}
    	if( aParam.has( "equivalentL" ) )
    	{
    		fEquivalentL = aParam["equivalentL"]().as_double();
    		fGeneratingTF = true;
    	}
    	if( aParam.has( "equivalentC" ) )
    	{
    		fEquivalentC = aParam["equivalentC"]().as_double();
    		fGeneratingTF = true;
    	}
    	if( aParam.has( "TFBins" ) )
    	{
    		fTFBins = aParam["TFBins"]().as_int();
    	}
    	if( aParam.has( "FreqRangeCenter" ) )
    	{
    		fFreqRangeCenter = aParam["FreqRangeCenter"]().as_double();
    	}

    	if (fGeneratingTF) GenerateTransferFunction();
		
    	return true;
    }

    void EquivalentCircuit::GenerateTransferFunction()
    {
    	double freq_min = (1.-0.1)*fFreqRangeCenter;
    	double freq_max = (1.+0.1)*fFreqRangeCenter;
    	double PI = 3.1415926;
    	double f_res = 1./2./PI/sqrt(fEquivalentL*fEquivalentC);
    	double Q = 1./fEquivalentR*sqrt(fEquivalentL/fEquivalentC);
    	//printf("f_res and Q are %g and %g\n", f_res, Q);

    	double freq[fTFBins];
    	double V_re[fTFBins];
    	double V_im[fTFBins];
    	FILE * fTFout = fopen("output/Generated_TF.txt", "w");
    	//	fprintf(fTFout,"# RLC circuit analytic TF with f_res = 1.0 GHz, Q = 1000\n#        Freq [GHz]     vRe []     vIm []\n");
    	fprintf(fTFout,"#Transfer function: Frequency, magnitude, phase, vRe, vIm\n");

    	for (int counts=0; counts<fTFBins; counts++)
    	{
    		freq[counts] = freq_min + counts*(freq_max - freq_min)/(1.0 * fTFBins);

    		double Z_LC = (2*PI*freq[counts]*fEquivalentL - 1./(2*PI*freq[counts])/fEquivalentC);
    		double phase = atan(  Z_LC / fEquivalentR);
    		double VoutByVin = fEquivalentR/sqrt( fEquivalentR*fEquivalentR + Z_LC*Z_LC);
    		V_re[counts] = VoutByVin * cos(phase);
    		V_im[counts] = VoutByVin * sin(phase);
    		double f_GHz = freq[counts] * 1.0e-9;
//              fprintf(fTFout,"%16g%16.8g%16.8g\n", f_GHz,V_re[counts], V_im[counts]);
    		fprintf(fTFout,"%e,%e,%e,%e,%e\n", freq[counts], VoutByVin,phase,V_re[counts], V_im[counts]);
    		const std::complex<double> temp(V_re[counts], V_im[counts]);
    		tfArray.push_back(temp);
    	}
    	initialFreq = freq[0];
    	fclose(fTFout);
    }



} /* namespace locust */

