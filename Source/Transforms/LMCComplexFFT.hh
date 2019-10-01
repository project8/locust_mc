/*
 * LMCComplexFFT.hh
 *
 *  Created on: Sep 30, 2019
 *      Author: P. T. Surukuchi
 */

#ifndef LMCCOMPLEXFFT_HH_
#define LMCCOMPLEXFFT_HH_

#include <fftw3.h>
#include <math.h>
#include "param.hh"


namespace locust
{
    /*!
     @class ComplexFFT
     @author P. T. Surukuchi
     @brief Class to handle one dimenensional discrete Fourier transforms (DFT) and discrete Fourier transforms (IDFT) using FFTW
     @details
     Currently being used for tranforming transfer functions to FIR hence only the reverse FFT functionality has been implemented. If Forward FFT implemented, other transform classes like Hilbert transform and /Generators/LMC*PassFilterFFTGenerator could use the functionality from here.
     By default uses
     
     Available configuration options:
     -"Transform-flag": string -- The flags argument is FFTW_MEASURE (default) or FFTW_ESTIMATE
     -"use-wisdom": bool -- Option to use FFTW wisdom file if it already exists to improve FFT performance
     -"wisdom-filename": string -- The wisdom file name to use
     
     Information from FFTW website (http://www.fftw.org/fftw3_doc/Complex-One_002dDimensional-DFTs.html#Complex-One_002dDimensional-DFTs)
     FFTW_MEASURE instructs FFTW to run and measure the execution time of several FFTs in order to find the best way to compute the transform of size n. This process takes some time (usually a few seconds), depending on your machine and on the size of the transform.
     FFTW_ESTIMATE, on the contrary, does not run any computation and just builds a reasonable plan that is probably sub-optimal. In short, if your program performs many transforms of the same size and initialization time is not important, use FFTW_MEASURE; otherwise use the estimate.
     */
    class ComplexFFT
    {
        
    public:
        ComplexFFT();
        virtual ComplexFFT();
        bool Configure( const scarab::param_node& aNode);
        
        enum Domain
        {
            kNone,
            kFrequencyDomain,
            kTimeDomain
        };
        
    private:
        
        fftw_complex* ReverseTransform();
        fftw_complex* ForwardTransform();
    };
    
} /* namespace locust */

#endif /* LMCHILBERTTRANSFORM_HH_ */
