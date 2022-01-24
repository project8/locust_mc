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


    EquivalentCircuit::EquivalentCircuit() {}
    EquivalentCircuit::~EquivalentCircuit() {}

    void EquivalentCircuit::GenerateTransferFunction(double R, double L, double C, int nBins_config = 4000, double FreqRangeCenter = 1.0e9){
	double freq_min = (1.-0.1)*FreqRangeCenter;
	double freq_max = (1.+0.1)*FreqRangeCenter;
	double PI = 3.1415926;
	nbins = nBins_config;
	double f_res = 1./2./PI/sqrt(L*C);
	double Q = 1./R*sqrt(L/C);
	//printf("f_res and Q are %g and %g\n", f_res, Q);

	double freq[nbins];
	double V_re[nbins];
	double V_im[nbins];
	FILE * fTFout = fopen("Generated_TF.txt", "w");
	fprintf(fTFout,"#Transfer function: Frequency, magnitude, phase, vRe, vIm\n");
	for (int counts=0; counts<nbins; counts++){
		freq[counts] = freq_min + counts*(freq_max - freq_min)/(1.0 * nbins);

		double Z_LC = (2*PI*freq[counts]*L - 1./(2*PI*freq[counts])/C);
		double phase = atan(  Z_LC / R);
		double VoutByVin = R/sqrt( R*R + Z_LC*Z_LC);
		V_re[counts] = VoutByVin * cos(phase);
		V_im[counts] = VoutByVin * sin(phase);
		fprintf(fTFout,"%e,%e,%e,%e,%e\n", freq[counts], VoutByVin,phase,V_re[counts], V_im[counts]);
		const std::complex<double> temp(V_re[counts], V_im[counts]);
		tfArray.push_back(temp);
	}
	initialFreq = freq[0];
	fclose(fTFout);
    }



} /* namespace locust */

