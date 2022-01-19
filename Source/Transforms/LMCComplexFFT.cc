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
    fInputWindowArray(NULL),
    fOutputWindowArray(NULL),
    fInputSignalArray(NULL),
    fOutputSignalArray(NULL),
    fNShiftBins(2000),
    fForwardPlan(),
    fReversePlan()
    {
    }
    
    ComplexFFT::~ComplexFFT()
    {
        if (fInputWindowArray != NULL)
        {
            fftw_free(fInputWindowArray);
            fInputWindowArray = NULL;
        }
        if (fOutputWindowArray != NULL)
        {
            fftw_free(fOutputWindowArray);
            fOutputWindowArray = NULL;
        }
        if (fInputSignalArray != NULL)
        {
            fftw_free(fInputSignalArray);
            fInputSignalArray = NULL;
        }
        if (fOutputSignalArray != NULL)
        {
            fftw_free(fOutputSignalArray);
            fOutputSignalArray = NULL;
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


    bool ComplexFFT::SetupIFFTWindow(int size, double intialBinValue, double freqResolution, std::string windowname, double windowparam)
    {
        fSize=size;	
        fFreqResolution=freqResolution;
        fPreFilterBins=(int)(intialBinValue*1.0/freqResolution);
        fTotalWindowSize=fPreFilterBins+fZeroPaddingSize+fSize;
        fTimeResolution=1.0/(fTotalWindowSize*freqResolution);

        fInputWindowArray = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fTotalWindowSize);
        fOutputWindowArray = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fTotalWindowSize);

        fInputSignalArray = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fSize);
        fOutputSignalArray = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fSize);

        fWindowFunction.SetupWindow(windowname, windowparam);
        if(!fWindowFunction.GenerateWindowFunction(fTotalWindowSize, fPreFilterBins, fSize))
        {
            LERROR(lmclog,"Error generating window function");
            exit(-1);
        }
        return true;
    }

    bool ComplexFFT::SetupFFTWindow(int size, std::string windowname, double windowparam)
    {
        fSize=size; 
        fPreFilterBins=0;
        fTotalWindowSize=size;

        fInputWindowArray = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fTotalWindowSize);
        fOutputWindowArray = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fTotalWindowSize);

        fInputSignalArray = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fSize);
        fOutputSignalArray = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fSize);

        fWindowFunction.SetupWindow(windowname, windowparam);
        if(!fWindowFunction.GenerateWindowFunction(fTotalWindowSize, fPreFilterBins, fSize))
        {
            LERROR(lmclog,"Error generating window function");
            exit(-1);
        }
        return true;
    }

    bool ComplexFFT::GenerateFIR(int size, fftw_complex* in, fftw_complex* out)
    {
        if(!IsInitialized) return false;

        for (int i = 0; i < fTotalWindowSize; ++i)
        {
            if(i>fPreFilterBins && i<fSize+fPreFilterBins)
            {
                fInputWindowArray[i][0]=in[i-fPreFilterBins][0];
                fInputWindowArray[i][1]=in[i-fPreFilterBins][1];
            }
            else
            {   
                fInputWindowArray[i][0]=0.0;
                fInputWindowArray[i][1]=0.0;
            }
        }

        ReverseFFT(fTotalWindowSize, fInputWindowArray, fOutputWindowArray);
        
        if(!MakeFilterCausal(fOutputWindowArray))
        {
        	LERROR(lmclog,"Couldn't make FIR filter causal");
        	exit(-1);
        }

        for (int i = 0; i < fSize+2*fNShiftBins; ++i)
        {
            out[i][0]=fOutputWindowArray[i][0]*2/fTotalWindowSize;
        }
        
        return true;
    }

    bool ComplexFFT::ApplyWindowFunction(int size, fftw_complex* in, fftw_complex* out)
    {



        if(!fWindowFunction.IsWindowGenerated())
        {
            LERROR(lmclog,"Error applying window function, it has not been generated");
            exit(-1); 
        }

        for (int i = 0; i < size; ++i)
        {
            fInputSignalArray[i][0]=in[i][0];
            fInputSignalArray[i][1]=in[i][1];
        }

        double windowfactor = 0.0;
        for (int i = 0; i < size; ++i)
        { 
            windowfactor = fWindowFunction.GetWindowFunction()->at(i+fPreFilterBins);
            fOutputSignalArray[i][0]=windowfactor*fInputSignalArray[i][0];
            fOutputSignalArray[i][1]=windowfactor*fInputSignalArray[i][1];
        }

        for (int i = 0; i < size; ++i)
        {
            out[i][0]=fOutputSignalArray[i][0];
            out[i][1]=fOutputSignalArray[i][1];
        }

        return true;
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
