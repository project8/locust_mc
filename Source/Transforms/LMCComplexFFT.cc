/*
 * LMCComplexFFT.cc
 *
 *  Created on: Sep 30, 2019
 *      Author: P. T. Surukuchi
 */
#include <iostream>
#include <fstream>

#include "logger.hh"
#include "LMCComplexFFT.hh"
#include "LMCConst.hh"

namespace locust
{
    LOGGER( lmclog, "ComplexFFT" );

    ComplexFFT::ComplexFFT():
    IsInitialized(false),
    fTransformFlag("MEASURE"),
    fWisdomFilename("wisdom_complexfft.fftw3"),
    fSize(0),
    fTotalWindowSize(0),
    fZeroPaddingSize(100000),
    fWindowFunctionType(1),
    fInputArray(NULL),
    fOutputArray(NULL),
    fNShiftBins(2000),
    fForwardPlan(),
    fReversePlan()
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
        if(aParam.has("zero-padding-size"))
        {
            fZeroPaddingSize=aParam["zero-padding-size"]().as_int();
        }
        if( aParam.has( "shift-n-bins" ) )
        {
            fNShiftBins=aParam["shift-n-bins"]().as_int();
        }
        if(!fWindowFunction.Configure(aParam))
        {
            LERROR(lmclog,"Error configuring Window Function class");
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

    bool ComplexFFT::ForwardFFT(int size, fftw_complex* in, fftw_complex* out)
    {
        if(!IsInitialized) return false;
        
        fForwardPlan = fftw_plan_dft_1d(size, in, out, FFTW_FORWARD, fTransform);
        fftw_execute_dft(fForwardPlan, in, out);
        fftw_destroy_plan(fForwardPlan);
        return true;
    }
    

    bool ComplexFFT::ReverseFFT(int size, fftw_complex* in, fftw_complex* out)
    {
        if(!IsInitialized) return false;

        fReversePlan = fftw_plan_dft_1d(size, in, out, FFTW_BACKWARD, fTransform);
        fftw_execute_dft(fReversePlan, in, out);
        fftw_destroy_plan(fReversePlan);
        return true;
    }


    bool ComplexFFT::SetupIFFT(int size, double intialBinValue, double freqResolution)
    {
        fSize=size;
        fFreqResolution=freqResolution;
        fPreFilterBins=(int)(intialBinValue*1.0/freqResolution);
        fTotalWindowSize=fPreFilterBins+fZeroPaddingSize+fSize;
        fTimeResolution=1.0/(fTotalWindowSize*freqResolution);
        if(GenerateWindowFunction()==false)
        {
            LERROR(lmclog,"Error generating window function");
            exit(-1);
        }
        return true;
    }

    bool ComplexFFT::GenerateFIR(int size, fftw_complex* in, fftw_complex* out)
    {
        fInputArray = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fTotalWindowSize);
        fOutputArray = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fTotalWindowSize);
        if(!IsInitialized) return false;

        if(!ApplyWindowFunction(fInputArray))
        {
            LERROR(lmclog,"Couldn't apply windowing function");
            exit(-1);
        }

        ReverseFFT(fTotalWindowSize, fInputArray, fOutputArray);
        
        if(!MakeFilterCausal(fOutputArray))
        {
        	LERROR(lmclog,"Couldn't make FIR filter causal");
        	exit(-1);
        }
        for (int i = 0; i < 2*fSize+fNShiftBins; ++i)
        {
        }
        for (int i = 0; i < fSize+2*fNShiftBins; ++i)
        {
            out[i][0]=fOutputArray[i][0]*2/fTotalWindowSize;
        }
        return true;
    }

    bool ComplexFFT::GenerateWindowFunction()
    {
       return fWindowFunction.GenerateWindowFunction(fTotalWindowSize, fPreFilterBins, fSize);
    }
    
    bool ComplexFFT::ApplyWindowFunction(fftw_complex* in)
    {
        double windowfactor = 0.;

        for (int i = 0; i < fTotalWindowSize; ++i)
        {
            windowfactor = fWindowFunction.GetWindowFunction().at(i);
            fInputArray[i][0]=windowfactor*in[i][0];
            fInputArray[i][1]=windowfactor*in[i][1];
        }
    }

    bool ComplexFFT::MakeFilterCausal(fftw_complex* in)
    {
    std::vector<double> vectorToTransformReal;
    std::vector<double> vectorToTransformImag;
        for (int i = 0; i < fTotalWindowSize; ++i)
        {
        vectorToTransformReal.push_back(in[i][0]);
        vectorToTransformImag.push_back(in[i][1]);
    }
    std::rotate(vectorToTransformReal.rbegin(),vectorToTransformReal.rbegin()+fNShiftBins,vectorToTransformReal.rend());
    std::rotate(vectorToTransformImag.rbegin(),vectorToTransformImag.rbegin()+fNShiftBins,vectorToTransformImag.rend());
        for (int i = 0; i < fTotalWindowSize; ++i)
        {
        in[i][0]=vectorToTransformReal.at(i);
        in[i][1]=vectorToTransformReal.at(i);
    }
    return true;
    }

    double ComplexFFT::GetTimeResolution()
    {
        return fTimeResolution;
    }
    
    double ComplexFFT::GetFreqResolution()
    {
        return fFreqResolution;
    }

    int ComplexFFT::GetShiftNBins()
    {
	return fNShiftBins;
    }
} /* namespace locust */
