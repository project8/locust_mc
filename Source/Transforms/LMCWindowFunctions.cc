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
    fTotalWindowSize(0),
    fPreWindowBins(0),
    fSignalSize(0),
    fWindowFunctionType(1),
    fWindowFunction( 1 ),
    fTukeyWindowAlpha( 0.5 )
    {
    }
    
    WindowFunctions::~WindowFunctions()
    {
    }
    
    bool WindowFunctions::Configure(const scarab::param_node& aParam)
    {
        return true;
    }

    bool WindowFunctions::SetupWindow(std::string windowname, double windowparam)
    {
        if(windowname == "rectangular")
        {
            fWindowFunctionType = 0;
        }
        else if(windowname == "tukey")
        {
            fWindowFunctionType = 1;
            fTukeyWindowAlpha = windowparam;
            printf("In Setup function\n");
            printf("fWindowFunctionType is %d\n", fWindowFunctionType);
            printf("fTukeyWindowAlpha is %g\n", fTukeyWindowAlpha);
        }
        else if(windowname == "other")
        {
            fWindowFunctionType = 2;
        }
        else
        {
            return false;
        }
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

        ////////////////////////////
        //text file for testing.   
        std::ofstream windowfile;
        windowfile.open("windowfile.txt");
        for (int i = 0; i < fTotalWindowSize; ++i)
        {
            windowfile << fWindowFunction.at(i);
            windowfile << "\n";
        }
        windowfile.close();
        ////////////////////////////

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
        // Tukey window is defined as: (wikipedia)
        // w[n]=0.5*(1-cos((pi*2*n)/(alpha*N))); if 0<=n<alpha*N/2
        // w[n]= 1; alpha*N/2<=n<=N/2
        // w[N-n] = w[n] for 0 <= n <= N/2

        //double tukeyWindowAlpha = 0.5;
        int firstWindowFirstBin = fPreWindowBins;
        int midWindowFirstBin = fPreWindowBins+fTukeyWindowAlpha*fSignalSize/2.0;
        int midWindowFinalBin = fPreWindowBins+fSignalSize-fTukeyWindowAlpha*fSignalSize/2.0;
        int finalWindowFinalBin = fPreWindowBins+fSignalSize;
        int j = 0; // index for the window itself, inside the larger array            
        
        for (int i = 0; i < fTotalWindowSize; ++i)
        {
            if(i>fPreWindowBins && i<midWindowFirstBin)
            {
                fWindowFunction[i]=0.5*(1 - std::cos( (LMCConst::Pi()*2*j) / (fTukeyWindowAlpha*fSignalSize) ) );
                j++;
            }
            else if(i>=midWindowFirstBin && i<=midWindowFinalBin)
            {
                fWindowFunction[i]=1.0;
                j++;
            }
            else if(i>midWindowFinalBin && i<=finalWindowFinalBin)
            {
                fWindowFunction[i]=0.5*(1 - std::cos( (LMCConst::Pi()*2*(fSignalSize-j)) / (fTukeyWindowAlpha*fSignalSize) ) );
                j++;
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
