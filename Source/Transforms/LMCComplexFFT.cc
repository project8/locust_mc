/*
 * LMCComplexFFT.cc
 *
 *  Created on: Sep 30, 2019
 *      Author: P. T. Surukuchi
 */
#include <iostream>
#include <fstream>

#include "LMCComplexFFT.hh"
#include "LMCConst.hh"

namespace locust
{
    ComplexFFT::ComplexFFT():
    IsInitialized(false),
    fTransformFlag("MEASURE"),
    fWisdomFilename("wisdom_complexfft.fftw3"),
    fSize(0),
    fTotalWindowSize(0),
    fZeroPaddingSize(0),
    fWindowFunctionType(1),
    fWindowFunction(NULL),
    fInputArray(NULL),
    fOutputArray(NULL),
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
        fftw_destroy_plan(fReversePlan);
        fftw_destroy_plan(fForwardPlan);
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
        if(aParam.has("window-function-type"))
        {
            fWindowFunctionType=aParam["window-function-type"]().as_int();
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
    
    bool ComplexFFT::GenerateWindowFunction()
    {
	std::ofstream myfile;
	myfile.open ("WindowFunction.txt");
        if(fWindowFunctionType==1)
	{
	    //Tukey window is defined as 
	    //w[n]=0.5*(1+cos(pi(2n/(alpha*N)-1))); if 0<=n<alpha*N/2
	    //w[n]= 1; alpha*N/2<=n<=N(1-alpha/2)
	    //w[n]=0.5*(1+cos(pi(2n/(alpha*N)-2/alpha+1))) if N(1-alpha/2)<n<=N
	    double tukeyWindowAlpha = 0.5;
	    int firstWindowFirstBin = fPreFilterBins;
	    int midWindowFirstBin = fPreFilterBins+tukeyWindowAlpha*fSize/2.0;
	    int midWindowFinalBin = fPreFilterBins+fSize-tukeyWindowAlpha*fSize/2.0;
	    int finalWindowFinalBin = fPreFilterBins+fSize;
            for (int i = 0; i < fTotalWindowSize; ++i)
            {
	        fWindowFunction.push_back(0.0);
	    }

            for (int i = 0; i < fTotalWindowSize; ++i)
            {
		if(i>fPreFilterBins && i<midWindowFirstBin)
		{
	            fWindowFunction[i]=0.5*(1+std::cos(LMCConst::Pi()*(2*i/(tukeyWindowAlpha*fSize)-1)));
		}
		else if(i>=midWindowFirstBin && i<=midWindowFinalBin)
		{
		    fWindowFunction[i]=1.0;
		}
		else if(i>midWindowFinalBin && i<=finalWindowFinalBin)
		{
		    fWindowFunction[i]=0.5*(1+std::cos(LMCConst::Pi()*(2*i/(tukeyWindowAlpha*fSize)-1.0/tukeyWindowAlpha+1)));
		}
	        myfile<<i;
	        myfile<<",";
	        myfile<<fWindowFunction[i];
	        myfile<<"\n";
	    }
	}	
	myfile.close();
	return true;
    }

    bool ComplexFFT::ForwardFFT(int size, fftw_complex* in, fftw_complex* out)
    {
        if(!IsInitialized) return false;
        
        fInputArray = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size);
        fOutputArray = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size);
        
        fForwardPlan = fftw_plan_dft_1d(size,fInputArray,fOutputArray,FFTW_FORWARD,fTransform);
        fftw_execute_dft(fForwardPlan,in,out);
        return true;
    }
    
    bool ComplexFFT::SetupIFFT(int size, double intialBinValue, double resolution)
    {
        fSize=size;	
	fPreFilterBins=(int)(intialBinValue*1.0/resolution);
	fTotalWindowSize=fPreFilterBins+fZeroPaddingSize+fSize;
	if(GenerateWindowFunction()==false)
	{
	    exit(-1);
	}
	return true;
    }

    bool ComplexFFT::ReverseFFT(int size, fftw_complex* in, fftw_complex* out)
    {
	std::ofstream myfile;
	myfile.open ("example.txt");
        fInputArray = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fTotalWindowSize);
        fOutputArray = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fTotalWindowSize);
        if(!IsInitialized) return false;
        for (int i = 0; i < fTotalWindowSize; ++i)
        {
	    if(i>fPreFilterBins && i<fSize+fPreFilterBins)
	    {
              fInputArray[i][0]=fWindowFunction.at(i)*in[i-fPreFilterBins][0];
              fInputArray[i][1]=fWindowFunction.at(i)*in[i-fPreFilterBins][1];
	    }
	    else
	    {
	      fInputArray[i][0]=0.0;
	      fInputArray[i][1]=0.0;
	    }
        }
	fReversePlan= fftw_plan_dft_1d(fTotalWindowSize,fInputArray,fOutputArray,FFTW_BACKWARD,FFTW_ESTIMATE);
        
        fftw_execute(fReversePlan);
        for (int i = 0; i < fTotalWindowSize; ++i){
	   if(i<fSize) out[i][0]=fOutputArray[i][0]/fTotalWindowSize;
	   myfile<<i;
	   myfile<<",";
	   myfile<<fOutputArray[i][0]/fTotalWindowSize;
	   myfile<<"\n";
        }
	myfile.close();
        return true;
    }
    
} /* namespace locust */
