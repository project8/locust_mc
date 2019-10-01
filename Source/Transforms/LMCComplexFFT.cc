/*
 * LMCComplexFFT.cc
 *
 *  Created on: Sep 30, 2019
 *      Author: P. T. Surukuchi
 */

#include "LMCComplexFFT.hh"

namespace locust
{
    ComplexFFT::ComplexFFT():
    IsInitialized(false),
    fTransformFlag("MEASURE"),
    fWisdomFilename("wisdom_complexfft.fftw3")
    {
    }
    
    ComplexFFT::~ComplexFFT()
    {
        if (fInputArray != NULL)
        {
            fftw_free(fInputArray);
            fInputArray = NULL;
        }
        
        if (fOutputArray != NULL)
        {
            fftw_free(fOutputArray);
            fOutputArray = NULL;
        }
        fftw_destroy_plan(ReversePlan);
        fftw_destroy_plan(ForwardPlan);
    }
    
    bool ComplexFFT::Configure(const scarab::param_node& aParam)
    {
        if( aParam.has("transform-flag"))
        {
            fTransformFlag=aParam["transform-flag"]().as_bool();
        }
        if(aParam.has("use-wisdom"))
        {
            fUseWisdom=aParam["use-wisdom"]().as_bool();
        }
        if(aParam.has("wisdom-filename"))
        {
            fWisdomFilename=aParam["wisdom-filename"]().as_string();
        }
        
        if(fTransformFlag.compare("MEASURE"))
        {
            fTransform=Transform::measure;
        }
        else if(fTransformFlag.compare("ESTIMATE"))
        {
            fTransform=Transform::estimate;
        }
        else
        {
            return false;
        }
        IsInitialized=true;
        return true;
    }
    
    bool ComplexFFT::ForwardFFT()
    {
        if(IsInitialized) return false;
        
        fInputArray = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fSize);
        fOutputArray = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fSize);
        
        ForwardPlan = fftw_plan_dft_1d(fSize,fInputArray,fOutputArray,FFTW_FORWARD,fTransform);
        return true;
    }
    
    bool ComplexFFT::ReverseFFT()
    {
        if(IsInitialized) return false;
        
        fInputArray = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fSize);
        fOutputArray = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fSize);
        
        ReversePlan = fftw_plan_dft_1d(fSize,fInputArray,fOutputArray,FFTW_BACKWARD,fTransform);
        return true;
    }

} /* namespace locust */
