/*
 * LMCHilbertTransform.cc
 *
 *  Created on: May 20, 2019
 *      Author: pslocum
 */

#include "LMCHilbertTransform.hh"
#include "logger.hh"

namespace locust
{

	LOGGER( lmclog, "HilbertTransform" );


    HilbertTransform::HilbertTransform():
        fbufferMargin( 25 ),
		fbufferSize( 50 )
    {
    }

    HilbertTransform::~HilbertTransform()
    {
    }

    bool HilbertTransform::Configure(const scarab::param_node& aParam)
    {
        if(!fComplexFFT.Configure(aParam))
        {
            LERROR(lmclog,"Error configuring ComplexFFT class");
        }
    	if( aParam.has( "hilbert-buffer-margin" ) )
        {
    		fbufferMargin=aParam["hilbert-buffer-margin"]().as_int();
        }
       	if( aParam.has( "hilbert-buffer-size" ) )
        {
       		fbufferSize=aParam["hilbert-buffer-size"]().as_int();
        }
       	if (2*fbufferMargin > fbufferSize)
       	{
       		return false;
       	}
       	else
       	{
       		return true;
       	}
    }

    void HilbertTransform::SetBufferSize( int abufferSize )
    {
    	fbufferSize = abufferSize;
    }
    void HilbertTransform::SetBufferMargin( int abufferMargin )
    {
    	fbufferMargin = abufferMargin;
    }
    int HilbertTransform::GetBufferSize()
    {
    	return fbufferSize;
    }
    int HilbertTransform::GetBufferMargin()
    {
    	return fbufferMargin;
    }


    double* HilbertTransform::GetMagPhaseMean(std::deque<double> FieldBuffer, std::deque<double> FrequencyBuffer)
    {

    	fftw_complex* transformeddata = Transform( FieldBuffer );
    	unsigned hilbertindex = fbufferMargin;
    	double* magphasemean = new double[3];

    	double mag = pow(transformeddata[hilbertindex][0]*transformeddata[hilbertindex][0] + transformeddata[hilbertindex][1]*transformeddata[hilbertindex][1], 0.5);
    	double mean = 0.;  // this is set to zero in HilbertTransform::Transform().
    	double phase = GetPhase(transformeddata[hilbertindex][0], transformeddata[hilbertindex][1], mean);

    	magphasemean[0] = mag;
    	magphasemean[1] = phase;
    	magphasemean[2] = mean;

    	delete[] transformeddata;

    	return magphasemean;
    }


    double HilbertTransform::GetPhase( double VI, double VQ, double VMean)
    {

    	double phase = 0.;
    	if (fabs(VI) > 0.) phase = atan(VQ/VI);

    	// atan(Q/I) only outputs 2 quadrants.  need all 4.
    	phase += QuadrantCorrection( VI, phase, VMean);

    	return phase;
    }


    double HilbertTransform::QuadrantCorrection( double VI, double HilbertPhase, double HilbertMean)
    {
    	double phasecorrection = 0.;
    	if (((HilbertPhase < 0.)&&(VI < HilbertMean)) || ((HilbertPhase > 0.)&&(VI < HilbertMean)))
    	{
    		phasecorrection = LMCConst::Pi();  // check IQ quadrant
    	}
    	return phasecorrection;
    }


    double HilbertTransform::GetMean( std::deque<double> FieldBuffer )
    {
    	double mean = 0.;
    	for (std::deque<double>::iterator it = FieldBuffer.begin(); it!=FieldBuffer.end(); ++it)
        {
    		mean += *it/(double)FieldBuffer.size();
        }
    	return mean;
    }


    double HilbertTransform::GetMean( fftw_complex* array, int IQ, int size )
    {
    	double mean = 0.;
    	for (unsigned i=0; i<size; i++)
        {
    		mean += array[i][IQ]/size;
        }
    	return mean;
    }


    double* HilbertTransform::GetSpan( fftw_complex* array, int IQ, int size )
    {
    	double* span = new double[2];
    	double max = -99.;
    	double min = 99.;
    	for (unsigned i=size/4; i<3*size/4; i++)
    	{
    		if (array[i][IQ] > max) max=array[i][IQ];
    		if (array[i][IQ] < min) min=array[i][IQ];
        }
    	span[0] = min;
    	span[1] = max;
    	return span;
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

        int i=0;
        for (std::deque<double>::iterator it = FieldBuffer.begin(); it!=FieldBuffer.end(); ++it)
        {
        	SignalComplex[i][0] = *it;
        	originaldata[i][0] = *it;
        	SignalComplex[i][1] = 0.;
        	originaldata[i][1] = 0.;
        	i += 1;
        }



        fComplexFFT.ForwardFFT(windowsize, SignalComplex, FFTComplex);


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


        fComplexFFT.ReverseFFT(windowsize, hilbert, SignalComplex);


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


        double* spanI = GetSpan(originaldata, 0, FieldBuffer.size());
        double* spanQ = GetSpan(originaldata, 1, FieldBuffer.size());

        for (int i = 0; i < windowsize; ++i)
        {
        	originaldata[i][0] -= (spanI[1]+spanI[0])/2.;
        	originaldata[i][1] -= (spanQ[1]+spanQ[0])/2.;
        }


        // dump to text file for debugging.
/*    	FILE *fp = fopen("Hilbertresults.txt", "w");
    	for (int i=0; i<windowsize; i++)
    	{
    		fprintf(fp, "%g %g %g %g\n", originaldata[i][0], originaldata[i][1], pow(originaldata[i][0]*originaldata[i][0] + originaldata[i][1]*originaldata[i][1], 0.5), GetPhase(originaldata[i][0], originaldata[i][1], 0.));
    	}
    	fclose (fp);
    	getchar();  // Control-C to quit.
*/

        delete[] hilbert;
        delete[] SignalComplex;
        delete[] FFTComplex;

        return originaldata;
   
    }



} /* namespace locust */
