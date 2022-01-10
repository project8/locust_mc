/*
 * LMCWindowFunctions.cc
 *
 *  Created on: Jan 10, 2022
 *      Author: A. B. Telles
 */
#include <iostream>
#include <fstream>

#include "logger.hh"
#include "LMCWindowFunctions.hh"
#include "LMCConst.hh"

namespace locust
{
    LOGGER( lmclog, "WindowFunctions" );

    WindowFunctions::WindowFunctions():
    IsInitialized(false),
    fTotalWindowSize(0),
    fPreWindowBins(0),
    fSignalSize(0),
    fWindowFunctionType(1),
    fWindowFunction( 1 )
    {
    }
    
    WindowFunctions::~WindowFunctions()
    {
    }
    
    bool WindowFunctions::Configure(const scarab::param_node& aParam)
    {
        if(aParam.has("window-function-type"))
        {
            if(aParam["window-function-type"]().as_string() == "rectangular")
            {
                fWindowFunctionType = 0;
            }
            if(aParam["window-function-type"]().as_string() == "tukey")
            {
                fWindowFunctionType = 1;
            }
            if(aParam["window-function-type"]().as_string() == "other")
            {
                fWindowFunctionType = 2;
            }
        }

        IsInitialized=true;
        return true;
    }
    
    bool WindowFunctions::GenerateWindowFunction(int windowsize, int prewindowbins, int signalsize)
    {
        fTotalWindowSize = windowsize;
        fPreWindowBins = prewindowbins;
        fSignalSize = signalsize;

        for (int i = 0; i < fTotalWindowSize; ++i)
        {
            fWindowFunction.push_back(0.0);
        }

        if(fWindowFunctionType==0)
        {
            GenerateRectangularWindow();
	    }
        else if(fWindowFunctionType==1)
        {
            GenerateTukeyWindow();
        }
        else if(fWindowFunctionType==2)
        {
            GenerateOtherWindow();
        }	
        else
        {
            return false;
        }
        return true;
    }

    void WindowFunctions::GenerateRectangularWindow()
    {
        for (int i = 0; i < fTotalWindowSize; ++i)
        {
            if(i>fPreWindowBins && i<=fPreWindowBins+fSignalSize)
            {
                fWindowFunction[i]=1.0;
            }
            else
            {
                fWindowFunction[i]=0.0;
            }
        }
    }

    void WindowFunctions::GenerateTukeyWindow()
    {
        //Tukey window is defined as 
        //w[n]=0.5*(1+cos(pi(2n/(alpha*N)-1))); if 0<=n<alpha*N/2
        //w[n]= 1; alpha*N/2<=n<=N(1-alpha/2)
        //w[n]=0.5*(1+cos(pi(2n/(alpha*N)-2/alpha+1))) if N(1-alpha/2)<n<=N

        double tukeyWindowAlpha = 0.5;
        int firstWindowFirstBin = fPreWindowBins;
        int midWindowFirstBin = fPreWindowBins+tukeyWindowAlpha*fSignalSize/2.0;
        int midWindowFinalBin = fPreWindowBins+fSignalSize-tukeyWindowAlpha*fSignalSize/2.0;
        int finalWindowFinalBin = fPreWindowBins+fSignalSize;
        
        for (int i = 0; i < fTotalWindowSize; ++i)
        {
           fWindowFunction.push_back(0.0);
        }
        for (int i = 0; i < fTotalWindowSize; ++i)
        {
            if(i>firstWindowFirstBin && i<midWindowFirstBin)
            {
                fWindowFunction[i]=0.5*(1+std::cos(LMCConst::Pi()*(2*i/(tukeyWindowAlpha*fSignalSize)-1)));
            }
            else if(i>=midWindowFirstBin && i<=midWindowFinalBin)
            {
                fWindowFunction[i]=1.0;
            }
            else if(i>midWindowFinalBin && i<=finalWindowFinalBin)
            {
                fWindowFunction[i]=0.5*(1+std::cos(LMCConst::Pi()*(2*i/(tukeyWindowAlpha*fSignalSize)-1.0/tukeyWindowAlpha+1)));
            }
        }
    }

    void WindowFunctions::GenerateOtherWindow()
    {
        for (int i = 0; i < fTotalWindowSize; ++i)
        {
            fWindowFunction[i]=0.0;
        }
    }

    std::vector<double> WindowFunctions::GetWindowFunction()
    {
        return fWindowFunction;
    }
    
} /* namespace locust */
