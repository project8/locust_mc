/*
 * LMCHilbertTransform.cc
 *
 *  Created on: May 20, 2019
 *      Author: pslocum
 */

#include "LMCHilbertTransform.hh"


namespace locust
{

    HilbertTransform::HilbertTransform()
    {
    }

    HilbertTransform::~HilbertTransform()
    {
    }



/*
    double* HilbertTransform::GetMagPhaseMean(std::deque<double> FieldBuffer, std::deque<double> FrequencyBuffer, unsigned edge_margin, double AcquisitionRate)
    {

    fftw_complex* transformeddata = Transform( FieldBuffer );
    double* frequencydata = GetFrequencyData( FrequencyBuffer );
    unsigned hilbertindex = edge_margin;
    double* magphasemean = new double[3];

    double mean = GetMean(FieldBuffer);
    double mag = pow(transformeddata[hilbertindex][0]*transformeddata[hilbertindex][0] + transformeddata[hilbertindex][1]*transformeddata[hilbertindex][1], 0.5);
    double phase = GetPhase(transformeddata[hilbertindex][0], transformeddata[hilbertindex][1], mean);

    magphasemean[0] = mag;
    magphasemean[1] = phase;
    magphasemean[2] = mean;

    // extrapolate to edge of window so we know the field phase on arrival at patch.
    for (unsigned i=hilbertindex; i>0; i--)
      magphasemean[1] -= 2.*LMCConst::Pi()*frequencydata[i]/AcquisitionRate;

    delete[] transformeddata;
    delete[] frequencydata;

    return magphasemean;
    }
*/



    double* HilbertTransform::GetMagPhaseMean(std::deque<double> FieldBuffer, std::deque<double> FrequencyBuffer, unsigned edge_margin, double AcquisitionRate)
    {

    fftw_complex* transformeddata = Transform( FieldBuffer );
    double* frequencydata = GetFrequencyData( FrequencyBuffer );
    unsigned hilbertindex = edge_margin;
    double* magphasemean = new double[3];

    double mag = pow(transformeddata[hilbertindex][0]*transformeddata[hilbertindex][0] + transformeddata[hilbertindex][1]*transformeddata[hilbertindex][1], 0.5);
    double mean = GetMean( FieldBuffer );
    double phase = GetPhase(transformeddata[hilbertindex][0], transformeddata[hilbertindex][1], mean);

    magphasemean[0] = mag;
    magphasemean[1] = phase;
    magphasemean[2] = mean;

    // extrapolate to edge of window so we know the field phase on arrival at patch.
    for (unsigned i=hilbertindex; i>0; i--)
      magphasemean[1] -= 2.*LMCConst::Pi()*frequencydata[i]/AcquisitionRate;

    delete[] transformeddata;
    delete[] frequencydata;

    return magphasemean;
    }


    double HilbertTransform::GetPhase( double VI, double VQ, double VMean)
    {

    double phase = 0.;
    if (fabs(VI) > 0.)
      phase = atan(VQ/VI);

    // atan(Q/I) only outputs 2 quadrants.  need all 4.
    phase += QuadrantCorrection( VI, phase, VMean);

    return phase;
    }


    double HilbertTransform::QuadrantCorrection( double VI, double HilbertPhase, double HilbertMean)
    {
    double phasecorrection = 0.;
    if (((HilbertPhase < 0.)&&(VI < HilbertMean)) ||
        ((HilbertPhase > 0.)&&(VI < HilbertMean)))
        phasecorrection = LMCConst::Pi();  // check IQ quadrant
    else phasecorrection = 0.;
    return phasecorrection;
    }


    double HilbertTransform::GetMean( std::deque<double> FieldBuffer )
    {
    double mean = 0.;
    for (std::deque<double>::iterator it = FieldBuffer.begin(); it!=FieldBuffer.end(); ++it)
          {
          mean += *it/FieldBuffer.size();
          }
    return mean;
    }


    double* HilbertTransform::GetFrequencyData( std::deque<double> FrequencyBuffer )
    {

    // unpack buffer into array and return it.
    double* frequencydata = new double[FrequencyBuffer.size()];

        int i=0;
        for (std::deque<double>::iterator it = FrequencyBuffer.begin(); it!=FrequencyBuffer.end(); ++it)
          {
          frequencydata[i] = *it;
          i += 1;
          }

    return frequencydata;
    }


    fftw_complex* HilbertTransform::Transform(std::deque<double> FieldBuffer)
    {

        int windowsize=FieldBuffer.size();

        fftw_complex *originaldata;
        originaldata = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * windowsize);
        fftw_complex *SignalComplex;
        SignalComplex = (fftw_complex*)fftw_malloc( sizeof(fftw_complex) * windowsize );
        fftw_complex *FFTComplex;
        FFTComplex = (fftw_complex*)fftw_malloc( sizeof(fftw_complex) * windowsize );
        fftw_complex *hilbert;
        hilbert = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * windowsize);

        fftw_plan ForwardPlan;
        ForwardPlan = fftw_plan_dft_1d(windowsize, SignalComplex, FFTComplex, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_plan ReversePlan;
        ReversePlan = fftw_plan_dft_1d(windowsize, hilbert, SignalComplex, FFTW_BACKWARD, FFTW_ESTIMATE);



        int i=0;
        for (std::deque<double>::iterator it = FieldBuffer.begin(); it!=FieldBuffer.end(); ++it)
          {
          SignalComplex[i][0] = *it;
          originaldata[i][0] = *it;
          SignalComplex[i][1] = 0.;
          originaldata[i][1] = 0.;
          i += 1;
          }



        fftw_execute(ForwardPlan); // SignalComplex->FFTComplex


	// do the phase shifts
  	for (int i = 0; i < windowsize; ++i)
          {
      	  if (i>windowsize/2+1)  // negative frequencies
            {     // (I+iQ)(-i)(-1) = iI-Q
            hilbert[i][0] = -FFTComplex[i][1]; 
            hilbert[i][1] = FFTComplex[i][0];
            }
          else  // positive frequencies
            {     //  (I+iQ)(-i)(1) = -iI+Q
            hilbert[i][0] = FFTComplex[i][1];
            hilbert[i][1] = -FFTComplex[i][0]; 
            }
          }

      fftw_execute(ReversePlan); // hilbert->SignalComplex


      for (int i = 0; i < windowsize; i++)  // normalize with 1/N
        {
        SignalComplex[i][0] /= windowsize;
        SignalComplex[i][1] /= windowsize;
        }



        // construct transformed time series.
        //    	    Y(t) = y(t) + j h(t) = originaldata + j(hilbert0 + j hilbert1)
      for (int i = 0; i < windowsize; ++i)
        {
        originaldata[i][0] = originaldata[i][0] - SignalComplex[i][1];
        originaldata[i][1] = SignalComplex[i][0];
        }

      fftw_destroy_plan(ForwardPlan);
      fftw_destroy_plan(ReversePlan);
      delete[] hilbert;
      delete[] SignalComplex;
      delete[] FFTComplex;



    // dump to text file for debugging.
/*        
    FILE *fp = fopen("Hilbertresults.txt", "w");
    for (int i=0; i<windowsize; i++)
    {
          fprintf(fp, "%g %g %g\n", originaldata[i][0], originaldata[i][1], pow(originaldata[i][0]*originaldata[i][0] + originaldata[i][1]*originaldata[i][1], 0.5));
    }
    fclose (fp);
    getchar();  // Control-C to quit.
*/

    return originaldata;
   
    }



} /* namespace locust */
