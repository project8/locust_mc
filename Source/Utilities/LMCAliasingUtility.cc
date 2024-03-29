/*
 * LMCAliasingUtility.cc
 *
 *  Created on: Feb 10, 2023
 *      Author: pslocum
 */

#include "LMCAliasingUtility.hh"
#include "logger.hh"

using namespace scarab;

namespace locust
{
    LOGGER( testlog, "AliasingUtility" );

    AliasingUtility::AliasingUtility()
    {
    }

    AliasingUtility::~AliasingUtility()
    {
    }

	bool AliasingUtility::Configure()
	{
		return true;
	}

    bool AliasingUtility::CheckAliasing( double RF, double LO, double fs, double dr )
    {
    	// RF is the rf frequency (Hz)
    	// LO is the local oscillator frequency (Hz)
    	// fs is the sampling frequency (MHz)
    	// dr is the decimation rate, presently hard-wired to 10.0 (unitless factor).
    	// fsFast (Hz) is the fast sampling rate used to help the LPF work correctly.

    	double fsFast = fs * 1.e6 * dr;
    	double nyquistFrequency = fsFast / dr / 2.;

    	int sign = 1; // start with upper sideband.

    	bool bPass = true;

    	for (int j=0; j<2; j++)  // choose upper or lower mixing sideband
    	{
    	    if (j==1) sign = -1;
    	    for (int i=0; i<3; i++) // first 3 harmonics
    	    {
    	    double freq = i*RF - sign*LO;
    	    int nwin = (round)(freq/fsFast);

    	    double alias = sign*(freq - nwin*fsFast);
    	    if (j==0)
    	    {
    	    	printf("%5d*RF - LO = %10.4g and alias is %2d*(%10.4g - %4d*%g)=%10.4g Hz\n", i, freq, sign, freq, nwin, fsFast, alias);
    	    	if ((i!=1) && (fabs(alias) < 2.*nyquistFrequency))
    	    	{
    	    		LERROR( testlog, "Aliased frequency " << fabs(alias) << " is below 2X Nyquist frequency " << 2.*nyquistFrequency );
    	    		bPass = false;
    	    	}
    	    }
    	    if (j==1)
    	    {
    	    	printf("%5d*RF + LO = %10.4g and alias is %2d*(%10.4g - %4d*%g)=%10.4g Hz\n", i, freq, sign, freq, nwin, fsFast, alias);
    	    	if (fabs(alias) < 2.*nyquistFrequency)
    	    	{
    	    		LERROR( testlog, "Aliased frequency " << fabs(alias) << " is below 2X Nyquist frequency " << 2.*nyquistFrequency );
    	    		bPass = false;
    	    	}
    	    }

    	    } // i
    	} // j

        LPROG( testlog, "\nSummary:");
    	LPROG( testlog, "Presently hard-wired decimation rate is " << dr );
    	LPROG( testlog, "acquisition-rate is " << fs << " MHz" );
    	LPROG( testlog, "lo-frequency is " << LO );
    	LPROG( testlog, "rf-frequency is " << RF );

    	if (bPass)
    	{
    		LPROG(testlog,"1*RF-LO is in the measurement window; other harmonics are not in the window.");
    		LPROG(testlog,"Aliasing checks have passed.");
    	}
    	else
    	{
    		LWARN(testlog,"There appears to be some potential HF aliasing.")
    	}


    	return bPass;
    }

} /* namespace locust */
