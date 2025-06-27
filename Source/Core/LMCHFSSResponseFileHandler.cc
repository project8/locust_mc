/*
 * LMCHFSSFileHandler.cc
 *
 *  Created on: May 11, 2018
 *      Author: P. T. Surukuchi
 */

#include "LMCConst.hh"
#include "LMCException.hh"
#include "LMCHFSSResponseFileHandler.hh"

#include <stdio.h>
#include "logger.hh"

#include<unistd.h>
#include <iomanip> 

namespace locust
{
    LOGGER( lmclog, "HFSSResponseFileHandlerCore" );
    
    HFSSResponseFileHandlerCore::HFSSResponseFileHandlerCore():
    fCharacteristicImpedance(1.0),
    fNSkips(1),
    fNModes(2),
    fComplexFFT(),
    fHFSSFiletype(""),
    fIsFIRCreated(false),
    fWindowName("tukey"),
    fWindowParam(0.5),
    fPrintFIR ( false ),
    fConvertStoZ ( false ),
    fOutputPath( TOSTRING(PB_OUTPUT_DIR))
    {
    }
    
    HFSSResponseFileHandlerCore::~HFSSResponseFileHandlerCore()
    {
    }
    
    bool HFSSResponseFileHandlerCore::Configure(const scarab::param_node& aParam)
    {
        if( aParam.has( "print-fir-debug" ) )
    	{
            fPrintFIR = aParam["print-fir-debug"]().as_bool();
    	}
        if( aParam.has( "n-modes" ) )
        {
            fNModes = aParam["n-modes"]().as_int();
        }

        DimensionMultiMode( fNModes );

        if( aParam.has( "convert-sparams-to-z"))
        {
        	fConvertStoZ = aParam["convert-sparams-to-z"]().as_bool();
        }

        if( aParam.has( "characteristic-impedance"))
        {
            fCharacteristicImpedance = aParam["characteristic-impedance"]().as_double();
            LPROG( lmclog, "Characteristic impedance has been changed to " << fCharacteristicImpedance);
        }

    	if ( aParam.has( "output-path" ) )
    	{
    		fOutputPath = aParam["output-path"]().as_string();
    	}

        return true;
    }

    bool HFSSResponseFileHandlerCore::DimensionMultiMode( int nModes )
    {
        fFilterComplexArray.resize(2); // TE or TM
        fTFNBinsArray.resize(2);
        fFIRNBinsArray.resize(2);
        fResolutionArray.resize(2);
        fIsFIRCreatedArray.resize(2);
        fEfield.resize(2);
        fBfield.resize(2);
        for ( unsigned bTE=0; bTE<2; bTE++)
        {
            fFilterComplexArray[bTE].resize(nModes);
            fTFNBinsArray[bTE].resize(nModes);
            fFIRNBinsArray[bTE].resize(nModes);
            fResolutionArray[bTE].resize(nModes);
            fIsFIRCreatedArray[bTE].resize(nModes);
            fEfield[bTE].resize(nModes);
            fBfield[bTE].resize(nModes);
            for (unsigned l=0; l<nModes; l++)
            {
                fFilterComplexArray[bTE][l].resize(nModes);
                fTFNBinsArray[bTE][l].resize(nModes);
                fFIRNBinsArray[bTE][l].resize(nModes);
                fResolutionArray[bTE][l].resize(nModes);
                fIsFIRCreatedArray[bTE][l].resize(nModes);
                fEfield[bTE][l].resize(nModes);
                fBfield[bTE][l].resize(nModes);
                for (unsigned m=0; m<nModes; m++)
                {
                    fFilterComplexArray[bTE][l][m].resize(nModes);
                    fTFNBinsArray[bTE][l][m].resize(nModes);
                    fFIRNBinsArray[bTE][l][m].resize(nModes);
                    fResolutionArray[bTE][l][m].resize(nModes);
                    fIsFIRCreatedArray[bTE][l][m].resize(nModes);
                    fEfield[bTE][l][m].resize(nModes);
                    fBfield[bTE][l][m].resize(nModes);
                }
            }
        }

    	return true;
    }
    
    bool HFSSResponseFileHandlerCore::ends_with(const std::string &str, const std::string &suffix)
    {
        //copied from https://stackoverflow.com/a/20446239
        return str.size() >= suffix.size() &&
        str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
    }
    
    bool HFSSResponseFileHandlerCore::ReadHFSSFile(int bTE, int l, int m, int n)
    {
        return true;
    }

    // OPTION B (DO NOT DELETE)
    std::pair<double,double> HFSSResponseFileHandlerCore::ConvolveWithComplexFIRFilterArray(int bTE, int l, int m, int n, std::deque<double> inputBuffer, double t)
    {   
        std::array<double, 2> lastEfield = {0.0, 0.0};
        std::array<double, 2> lastBfield = {0.0, 0.0};

        double fBFactor = fFilterComplexArray[bTE][l][m][n][fFIRNBinsArray[bTE][l][m][n] - 1][0];
        double wprime = fFilterComplexArray[bTE][l][m][n][fFIRNBinsArray[bTE][l][m][n] - 1][1];
        double dt = fResolutionArray[bTE][l][m][n];

        if (fEfield[bTE][l][m][n].size() > 2)
        {
            lastEfield = fEfield[bTE][l][m][n].back();
            lastBfield = fBfield[bTE][l][m][n].back();
        }
        else
        {
            Jm0 = inputBuffer[0];
            LPROG("Jm0: " << std::setprecision(16) << Jm0);
        }
        // LPROG("lastEfield: " << lastEfield[0] << " + i*" << lastEfield[1]);
        // LPROG("lastBfield: " << lastBfield[0] << " + i*" << lastBfield[1]);

        if(fFIRNBinsArray[bTE][l][m][n]<=0)
        {   
            LERROR(lmclog,"Number of bins in the complex array filter should be positive");
        }   

        int inputBufferSize = inputBuffer.size(); 
        
        fEfield[bTE][l][m][n].resize(inputBufferSize,std::array<double, 2>{0.0, 0.0});
        fBfield[bTE][l][m][n].resize(inputBufferSize,std::array<double, 2>{0.0, 0.0});

        for (int i=0; i<inputBufferSize; i++)
        {
            // Initialize all fields except the last to zero
            if (i < inputBufferSize - 1)
            {
                fEfield[bTE][l][m][n][i] = {0.0, 0.0};
                fBfield[bTE][l][m][n][i] = {0.0, 0.0};
            }
            // First field, initialization of final field
            if (i==0)
            {
                fEfield[bTE][l][m][n][0][0] += lastEfield[0]*fFilterComplexArray[bTE][l][m][n][3*1 + 2][0] 
                                                        - lastEfield[1]*fFilterComplexArray[bTE][l][m][n][3*1 + 2][1];
                fEfield[bTE][l][m][n][0][1] += lastEfield[0]*fFilterComplexArray[bTE][l][m][n][3*1 + 2][1]
                                                        + lastEfield[1]*fFilterComplexArray[bTE][l][m][n][3*1 + 2][0];
                fBfield[bTE][l][m][n][0][0] += lastBfield[0]*fFilterComplexArray[bTE][l][m][n][3*1 + 2][0]
                                                        - lastBfield[1]*fFilterComplexArray[bTE][l][m][n][3*1 + 2][1];
                fBfield[bTE][l][m][n][0][1] += lastBfield[0]*fFilterComplexArray[bTE][l][m][n][3*1 + 2][1]
                                                        + lastBfield[1]*fFilterComplexArray[bTE][l][m][n][3*1 + 2][0];

                // Last field gets initialized here at the start
                fEfield[bTE][l][m][n][inputBufferSize-1][0] = lastEfield[0] * fFilterComplexArray[bTE][l][m][n][3*inputBufferSize + 2][0] 
                                                                - lastEfield[1] * fFilterComplexArray[bTE][l][m][n][3*inputBufferSize + 2][1];
                fEfield[bTE][l][m][n][inputBufferSize-1][1] = lastEfield[0] * fFilterComplexArray[bTE][l][m][n][3*inputBufferSize + 2][1]
                                                                + lastEfield[1] * fFilterComplexArray[bTE][l][m][n][3*inputBufferSize + 2][0];
                fBfield[bTE][l][m][n][inputBufferSize-1][0] = lastBfield[0] * fFilterComplexArray[bTE][l][m][n][3*inputBufferSize + 2][0]
                                                                - lastBfield[1] * fFilterComplexArray[bTE][l][m][n][3*inputBufferSize + 2][1];
                fBfield[bTE][l][m][n][inputBufferSize-1][1] = lastBfield[0] * fFilterComplexArray[bTE][l][m][n][3*inputBufferSize + 2][1]
                                                                + lastBfield[1] * fFilterComplexArray[bTE][l][m][n][3*inputBufferSize + 2][0];
            }
            // Middle fields
            else if (i<inputBufferSize-1)
            {
                // Time evolution additions
                fEfield[bTE][l][m][n][i][0] += fEfield[bTE][l][m][n][i-1][0]*fFilterComplexArray[bTE][l][m][n][3*1 + 2][0]
                                                        - fEfield[bTE][l][m][n][i-1][1]*fFilterComplexArray[bTE][l][m][n][3*1 + 2][1];
                fEfield[bTE][l][m][n][i][1] += fEfield[bTE][l][m][n][i-1][0]*fFilterComplexArray[bTE][l][m][n][3*1 + 2][1]
                                                        + fEfield[bTE][l][m][n][i-1][1]*fFilterComplexArray[bTE][l][m][n][3*1 + 2][0];
                fBfield[bTE][l][m][n][i][0] += fBfield[bTE][l][m][n][i-1][0]*fFilterComplexArray[bTE][l][m][n][3*1 + 2][0]
                                                        - fBfield[bTE][l][m][n][i-1][1]*fFilterComplexArray[bTE][l][m][n][3*1 + 2][1];
                fBfield[bTE][l][m][n][i][1] += fBfield[bTE][l][m][n][i-1][0]*fFilterComplexArray[bTE][l][m][n][3*1 + 2][1]
                                                        + fBfield[bTE][l][m][n][i-1][1]*fFilterComplexArray[bTE][l][m][n][3*1 + 2][0];
            }
            // Final field gets special treatment
            if (i < inputBufferSize - 1)
            {
                // Impulse additions to final field
                fEfield[bTE][l][m][n][inputBufferSize-1][0] += inputBuffer[i]*fFilterComplexArray[bTE][l][m][n][3*(inputBufferSize - 1 - i)][0];
                fEfield[bTE][l][m][n][inputBufferSize-1][1] += inputBuffer[i]*fFilterComplexArray[bTE][l][m][n][3*(inputBufferSize - 1 - i)][1];
                fBfield[bTE][l][m][n][inputBufferSize-1][0] += inputBuffer[i]*fFilterComplexArray[bTE][l][m][n][3*(inputBufferSize - 1 - i) + 1][0];
                fBfield[bTE][l][m][n][inputBufferSize-1][1] += inputBuffer[i]*fFilterComplexArray[bTE][l][m][n][3*(inputBufferSize - 1 - i) + 1][1];
            }
                // LPROG("Test: " << inputBuffer[i]*fFilterComplexArray[bTE][l][m][n][3*(inputBufferSize - 1- i) + 1][0]);

            // Impulse additions to all field
            fEfield[bTE][l][m][n][i][0] += inputBuffer[i]*fFilterComplexArray[bTE][l][m][n][0][0];
            fEfield[bTE][l][m][n][i][1] += inputBuffer[i]*fFilterComplexArray[bTE][l][m][n][0][1];
            fBfield[bTE][l][m][n][i][0] += inputBuffer[i]*fFilterComplexArray[bTE][l][m][n][1][0];
            fBfield[bTE][l][m][n][i][1] += inputBuffer[i]*fFilterComplexArray[bTE][l][m][n][1][1];

            // LPROG(lmclog, "Efield[" << bTE << "][" << l << "][" << m << "][" << n << "] = " << fEfield[bTE][l][m][n][i][0] << " + i*" << fEfield[bTE][l][m][n][i][1]);
            // LPROG(lmclog, "Bfield[" << bTE << "][" << l << "][" << m << "][" << n << "]" << "[" << i << "] = " << fBfield[bTE][l][m][n][i][0] << " + i*" << fBfield[bTE][l][m][n][i][1]);

            // Initial phase term
            fEfield[bTE][l][m][n][i][0] += Jm0 * exp( - fBFactor * (t + (i+1)*dt)) * sin( wprime * (t + (i+1)*dt)) / wprime / LMCConst::EpsNull();
            fEfield[bTE][l][m][n][i][1] += -Jm0 * exp( - fBFactor * (t + (i+1)*dt)) * cos( wprime * (t + (i+1)*dt)) / wprime / LMCConst::EpsNull();

            // LPROG(lmclog, "Efield[" << bTE << "][" << l << "][" << m << "][" << n << "] = " << fEfield[bTE][l][m][n][i][0] << " + i*" << fEfield[bTE][l][m][n][i][1]);
            // LPROG(lmclog, "Bfield[" << bTE << "][" << l << "][" << m << "][" << n << "]" << "[" << i << "] = " << std::setprecision(16) << fBfield[bTE][l][m][n][i][0] << " + i*" << fBfield[bTE][l][m][n][i][1]);
            // LPROG("JdotE" << i << ": " << inputBuffer[i]);

            // if (inputBuffer[i] * fEfield[bTE][l][m][n][i][0] > 1e-24)
            // {
            //     LPROG("time: " << t + i*dt);
            //     std::abort();
            // }
        }

        return std::make_pair(fEfield[bTE][l][m][n].back()[0], fBfield[bTE][l][m][n].back()[0]);
    }  

    // OPTION A (DO NOT DELETE)
    // std::pair<double,double> HFSSResponseFileHandlerCore::ConvolveWithComplexFIRFilterArray(int bTE, int l, int m, int n, std::deque<double> inputBuffer, double t)
    // {   
    //     std::array<double, 2> lastEfield = {0.0, 0.0};
    //     std::array<double, 2> lastBfield = {0.0, 0.0};

    //     double fBFactor = fFilterComplexArray[bTE][l][m][n][fFIRNBinsArray[bTE][l][m][n] - 1][0];
    //     double wprime = fFilterComplexArray[bTE][l][m][n][fFIRNBinsArray[bTE][l][m][n] - 1][1];
    //     double dt = fResolutionArray[bTE][l][m][n];

    //     if (fEfield[bTE][l][m][n].size() > 2)
    //     {
    //         lastEfield = fEfield[bTE][l][m][n].back();
    //         lastBfield = fBfield[bTE][l][m][n].back();
    //     }
    //     else
    //     {
    //         Jm0 = inputBuffer[0];
    //         LPROG("Jm0: " << Jm0);
    //     }
    //     // LPROG("lastEfield: " << lastEfield[0] << " + i*" << lastEfield[1]);

    //     if(fFIRNBinsArray[bTE][l][m][n]<=0)
    //     {   
    //         LERROR(lmclog,"Number of bins in the complex array filter should be positive");
    //     }   

    //     int inputBufferSize = inputBuffer.size(); 
        
    //     fEfield[bTE][l][m][n].resize(inputBufferSize,std::array<double, 2>{0.0, 0.0});
    //     fBfield[bTE][l][m][n].resize(inputBufferSize,std::array<double, 2>{0.0, 0.0});

    //     for (int i=0; i<inputBufferSize; i++)
    //     {
    //         fEfield[bTE][l][m][n][i] = {0.0, 0.0};
    //         fBfield[bTE][l][m][n][i] = {0.0, 0.0};

    //         for (int j=0; j<i+1; j++)
    //         {
    //             // Impulse additions to field
    //             fEfield[bTE][l][m][n][i][0] += inputBuffer[j]*fFilterComplexArray[bTE][l][m][n][3*(i-j)][0];
    //             fEfield[bTE][l][m][n][i][1] += inputBuffer[j]*fFilterComplexArray[bTE][l][m][n][3*(i-j)][1];
    //             fBfield[bTE][l][m][n][i][0] += inputBuffer[j]*fFilterComplexArray[bTE][l][m][n][3*(i-j)+1][0];
    //             fBfield[bTE][l][m][n][i][1] += inputBuffer[j]*fFilterComplexArray[bTE][l][m][n][3*(i-j)+1][1];
    //         }
    //         // LPROG(lmclog, "Efield[" << bTE << "][" << l << "][" << m << "][" << n << "] = " << fEfield[bTE][l][m][n][i][0] << " + i*" << fEfield[bTE][l][m][n][i][1]);
    //         // LPROG(lmclog, "Bfield[" << bTE << "][" << l << "][" << m << "][" << n << "] = " << fBfield[bTE][l][m][n][i][0] << " + i*" << fBfield[bTE][l][m][n][i][1]);

    //         // Time evolution additions to field
    //         fEfield[bTE][l][m][n][i][0] += lastEfield[0]*fFilterComplexArray[bTE][l][m][n][3*(i+1)+2][0]
    //                                                     - lastEfield[1]*fFilterComplexArray[bTE][l][m][n][3*(i+1)+2][1];
    //         fEfield[bTE][l][m][n][i][1] += lastEfield[0]*fFilterComplexArray[bTE][l][m][n][3*(i+1)+2][1]
    //                                                     + lastEfield[1]*fFilterComplexArray[bTE][l][m][n][3*(i+1)+2][0];
    //         fBfield[bTE][l][m][n][i][0] += lastBfield[0]*fFilterComplexArray[bTE][l][m][n][3*(i+1)+2][0]
    //                                                     - lastBfield[1]*fFilterComplexArray[bTE][l][m][n][3*(i+1)+2][1];
    //         fBfield[bTE][l][m][n][i][1] += lastBfield[0]*fFilterComplexArray[bTE][l][m][n][3*(i+1)+2][1]
    //                                                     + lastBfield[1]*fFilterComplexArray[bTE][l][m][n][3*(i+1)+2][0];

    //         // LPROG(lmclog, "Efield[" << bTE << "][" << l << "][" << m << "][" << n << "] = " << fEfield[bTE][l][m][n][i][0] << " + i*" << fEfield[bTE][l][m][n][i][1]);

    //         // Initial phase term
    //         fEfield[bTE][l][m][n][i][0] += Jm0 * exp( - fBFactor * (t + (i+1)*dt)) * sin( wprime * (t + (i+1)*dt)) / wprime / LMCConst::EpsNull();
    //         fEfield[bTE][l][m][n][i][1] += -Jm0 * exp( - fBFactor * (t + (i+1)*dt)) * cos( wprime * (t + (i+1)*dt)) / wprime / LMCConst::EpsNull();

    //         // LPROG(lmclog, "Efield[" << bTE << "][" << l << "][" << m << "][" << n << "] = " << fEfield[bTE][l][m][n][i][0] << " + i*" << fEfield[bTE][l][m][n][i][1]);
    //         // LPROG(lmclog, "Bfield[" << bTE << "][" << l << "][" << m << "][" << n << "] = " << fBfield[bTE][l][m][n][i][0] << " + i*" << fBfield[bTE][l][m][n][i][1]);
    //         // LPROG("JdotE" << i << ": " << inputBuffer[i]);

    //         // if (inputBuffer[i] * fEfield[bTE][l][m][n][i][0] > 1e-24)
    //         // {
    //         //     LPROG("time: " << t + i*dt);
    //         //     std::abort();
    //         // }
    //     }
    //     return std::make_pair(fEfield[bTE][l][m][n].back()[0], fBfield[bTE][l][m][n].back()[0]);
    // }  
    

    double HFSSResponseFileHandlerCore::QuadrantCorrection( double aRealValue, double aPhase)
    {
    	double phasecorrection = 0.;
    	if (aRealValue < 0)
    	{
    		phasecorrection = LMCConst::Pi();  // check IQ quadrant
    	}
    	return phasecorrection;
    }


    TFFileHandlerCore::TFFileHandlerCore():HFSSResponseFileHandlerCore(),
    fTFComplex(NULL),
    fFIRComplex(NULL),
    fInitialTFIndex(0.0),
    fTFBinWidth(100e6)
    {
    }
    
    TFFileHandlerCore::~TFFileHandlerCore()
    {
        if (fTFComplex != NULL)
        {
            fftw_free(fTFComplex);
            fTFComplex = NULL;
        }
        
        if (fFIRComplex != NULL)
        {
            fftw_free(fFIRComplex);
            fFIRComplex = NULL;
        }


        for( int bTE=0; bTE<2; bTE++)
        {
            for(int l=0; l<fNModes; l++)
            {
                for(int m=0; m<fNModes; m++)
                {
                    for(int n=0; n<fNModes; n++)
                    {
                    	fftw_free(fFilterComplexArray[bTE][l][m][n]);
                    	fFilterComplexArray[bTE][l][m][n] = NULL;
                    }
                }
            }
        }
    }

    bool TFFileHandlerCore::Configure(const scarab::param_node& aParam)
    {
        return true;
    }



    bool TFFileHandlerCore::ConvertTFtoFIR(int bTE, int l, int m, int n, std::vector<std::complex<double>> &tfArray)
    {
        if(fTFNBinsArray[bTE][l][m][n]<=0)
        {
            LERROR(lmclog,"The size of transfer function has to be positive integer");
            return false;
        }

        //Might need to be moved to a different function
        fTFComplex=(fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fTFNBinsArray[bTE][l][m][n]);

        for (int i = 0; i < fTFNBinsArray[bTE][l][m][n]; ++i)
        {
            fTFComplex[i][0]=tfArray.at(i).real();
            fTFComplex[i][1]=tfArray.at(i).imag();
        }

        // this length is somewhat arbitrary, but works for now
        // future improvement: make it adjust depending on length of FIR
        fFIRNBinsArray[bTE][l][m][n]=fTFNBinsArray[bTE][l][m][n]+ 2*fComplexFFT.GetShiftNBins();
        // Use TF bin width as given in .json config file
        fComplexFFT.SetupIFFTWindow(fTFNBinsArray[bTE][l][m][n],fInitialTFIndex,fTFBinWidth, fWindowName, fWindowParam);

        fFIRComplex=(fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fFIRNBinsArray[bTE][l][m][n]);
        fFilterComplexArray[bTE][l][m][n]=(fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fFIRNBinsArray[bTE][l][m][n]);
        fComplexFFT.ApplyWindowFunction(fTFNBinsArray[bTE][l][m][n], fTFComplex);
        fComplexFFT.GenerateFIR(fTFNBinsArray[bTE][l][m][n],fTFComplex,fFIRComplex);
        fResolutionArray[bTE][l][m][n]=fComplexFFT.GetTimeResolution();

        for (int i = 0; i < fFIRNBinsArray[bTE][l][m][n]; i++)
        {
            fFilterComplexArray[bTE][l][m][n][i][0] = fFIRComplex[i][0];
            fFilterComplexArray[bTE][l][m][n][i][1] = fFIRComplex[i][1];
        }

        fIsFIRCreatedArray[bTE][l][m][n]=true;

        LDEBUG( lmclog, "Finished IFFT to convert transfer function to FIR");
        return true;
    }
    
    bool TFFileHandlerCore::ReadHFSSFile(int bTE, int l, int m, int n)
    {

    	if(fIsFIRCreated)
    	{
    		return true;
    	}
        int tTFNBins = 0;
        if(!ends_with(fHFSSFilename,".txt"))
        {
            LERROR(lmclog,"The TF file " << fHFSSFilename.c_str() <<"doesn't end in .txt");
            return false;
        }
        double tfIndex;
        double tfRealValue;
        double tfImaginaryValue;
        std::vector<std::complex<double>> tfArray;
        std::fstream tfFile(fHFSSFilename.c_str(),std::ios::in);
        if (tfFile.fail())
        {
            LERROR(lmclog,"The TF file " << fHFSSFilename.c_str() <<" doesn't exist");
            return false;
        }
        int totalcount=0;
        
        while (!tfFile.eof()){
            std::string lineContent;
            while(std::getline(tfFile,lineContent))
            {
                if (lineContent.find('#')==std::string::npos)
                {
                    if (totalcount%fNSkips!=0){
                        ++totalcount;
                        continue;
                    }
                    std::string token;
                    std::stringstream ss(lineContent);
                    int wordCount=0;
                    while (ss >> token)
                    {
                        if(wordCount==0)tfIndex=std::stod(token);
                        else if(wordCount==1)tfRealValue=std::stod(token);
                        else if(wordCount==2)tfImaginaryValue=std::stod(token);
                        else
                        {
                            LERROR(lmclog, "There are more column than expected in the input TF file");
                            return false;
                        }
                        ++wordCount;
                    }
                    // The TF values from HFSS are in GHz, so need to convert to Hz
                    if(tTFNBins == 0)fInitialTFIndex=tfIndex*pow(10.0,9);
                    const std::complex<double> temp(tfRealValue,tfImaginaryValue);
                    tfArray.push_back(temp);
                    ++tTFNBins;
                }
            }
        }
        tfFile.close();
        LDEBUG( lmclog, "Finished reading transfer function file");

        ConvertAnalyticTFtoFIR(bTE, l, m, n, fInitialTFIndex, tfArray);

        if (fPrintFIR)
        {
            std::string modeIndexStr = std::to_string(bTE) + std::to_string(l) + std::to_string(m) + std::to_string(n);
            std::string fileName = fOutputPath + "/FIRhisto" + modeIndexStr + ".root";
            PrintFIR( fFilterComplexArray[bTE][l][m][n], fFIRNBinsArray[bTE][l][m][n], fileName );
            PrintFIR( fTFComplex, tTFNBins, fOutputPath+"/TFhisto.root");
            LPROG( lmclog, "Finished writing histos to output/FIRhisto.root and output/TFhisto.root");
            LPROG( lmclog, "Press Return to continue, or Cntrl-C to quit.");
            getchar();
        }

        return true;

    }


    bool TFFileHandlerCore::ConvertStoZ(std::vector<std::complex<double>> &tfArray, bool bConvert)
    {
    	if (bConvert)
    	{
    	    int count = tfArray.size();
    	    double maxZ = 0.;
    	    double maxZFrequency = 0.;
    	    double qInferred = -99.;
    	    double tFrequency = fInitialTFIndex;

    	    for (int i=0; i<count; i++)
    	    {
    	        double R = tfArray.at(i).real();
    	        double X = tfArray.at(i).imag();
    	        double tfI = (1.-R*R-X*X)/((1.-R)*(1.-R)+X*X);
    	        double tfQ = -(2.*X)/((1.-R)*(1.-R)+X*X);
    	        double tfMagSquared = tfI*tfI + tfQ*tfQ;
    	        if ( (tfMagSquared > maxZ) && (tFrequency > fInitialTFIndex) )
    	        {
    	            maxZFrequency = tFrequency;
    	            maxZ = tfMagSquared;
    	            qInferred = 0.;
    	        }
    	        else if ((tfMagSquared < 0.5*maxZ) && (qInferred == 0.))
    	        {
    	            double dFrequency = tFrequency - maxZFrequency;
    	            qInferred = maxZFrequency /  (2.* dFrequency);
    	        }

    	        tfArray[i] = std::complex<double>(tfI, tfQ);
    	        tFrequency += fTFBinWidth*1.e9;
    	    }
    	    for (int i=0; i<count; i++)
    	    {
    		    tfArray[i] /= pow(maxZ,0.5);
    		    tfArray[i] *= fCharacteristicImpedance;
    	    }
    	    LPROG(lmclog, "Resonant frequency is " << maxZFrequency << " and inferred Q is " << qInferred);
    	}

    	return true;
    }

    bool TFFileHandlerCore::ConvertAnalyticGFtoFIR( std::vector<std::vector<int>> aModeSet, std::vector<std::vector<std::pair<double,std::pair<double,double>>>> gfArray)
    {
        for (int mu=0; mu < aModeSet.size(); mu++)
        {
            bool bTE = aModeSet[mu][0];
            int l = aModeSet[mu][1];
            int m = aModeSet[mu][2];
            int n = aModeSet[mu][3];
            fFIRNBinsArray[bTE][l][m][n] = gfArray[mu].size();
            fResolutionArray[bTE][l][m][n] = gfArray[mu][0].first;
            fFilterComplexArray[bTE][l][m][n]=(fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fFIRNBinsArray[bTE][l][m][n]);


            for (int i = 0; i < fFIRNBinsArray[bTE][l][m][n]; i++)
            {
                fFilterComplexArray[bTE][l][m][n][i][0] = gfArray[mu][i].second.first;
                fFilterComplexArray[bTE][l][m][n][i][1] = gfArray[mu][i].second.second;
            }

            if (fPrintFIR)
            {
                std::string modeIndexStr = std::to_string(bTE) + std::to_string(l) + std::to_string(m) + std::to_string(n);
                std::string fileName = fOutputPath + "/FIRhisto" + modeIndexStr + ".root";
                PrintFIR( fFilterComplexArray[bTE][l][m][n], fFIRNBinsArray[bTE][l][m][n], fileName );
                LPROG( lmclog, "Finished writing histos to output/FIRhisto["<<bTE<<"]["<<l<<"]["<<m<<"]["<<n<<"].root");
                LPROG( lmclog, "Press Return to continue, or Cntrl-C to quit.");
                getchar();
            }
            fIsFIRCreatedArray[bTE][l][m][n]=true;
            LDEBUG( lmclog, "Finished populating FIR filter with Green's function.");
        }

        return true;
    }


    bool TFFileHandlerCore::ConvertAnalyticTFtoFIR(int bTE, int l, int m, int n, double initialFreq, std::vector<std::complex<double>> tfArray)
    {

	    //Replaces ReadHFSSFile() in the case where a Transfer Funciton has been generated analytically

        if(fIsFIRCreatedArray[bTE][l][m][n])
        {
            return true;
        }
        fTFNBinsArray[bTE][l][m][n]=0;

	    fTFNBinsArray[bTE][l][m][n]=tfArray.size();
	    fInitialTFIndex = initialFreq;
	    if(!ConvertTFtoFIR(bTE, l, m, n, tfArray))
	    { //bool determines if TF was generated dynamically
            return false;
        }
        fIsFIRCreatedArray[bTE][l][m][n]=true;
        return true;
    }    

	bool HFSSResponseFileHandlerCore::WriteRootHisto( std::vector<double> aFilter, int nBins, bool bIQ )
	{
#ifdef ROOT_FOUND
	    char fbuffer[60];
	    if (!bIQ)
	    {
	    	int a = sprintf(fbuffer, "Real");
	    }
	    else
	    {
	    	int a = sprintf(fbuffer, "Imag");
	    }
		fRootHistoWriter->OpenFile("UPDATE");
		const char *hName = fbuffer;
		TH1D* aHisto = new TH1D(hName, "Coefficients; index; coefficient", nBins, 0., (double)nBins);
		aHisto->SetDirectory(0);

		for (unsigned i=0; i<nBins; i++)
		{
			aHisto->SetBinContent(i+1, aFilter[i]);
		}

		fRootHistoWriter->Write1DHisto(aHisto);
		fRootHistoWriter->CloseFile();
		delete aHisto;
#endif
		return true;
	}

    
    void HFSSResponseFileHandlerCore::PrintFIR( std::vector<double> aFilter, int nBins, std::string filename )
    {
    	LDEBUG( lmclog, "Printing coefficients to file ... ");
    	std::string textFile(filename);
    	std::string suffix(".txt");
    	textFile.replace(textFile.find(".root"),5,suffix);
    	FILE * fFIRout = fopen(textFile.c_str(), "w");
    	for (int i = 0; i < nBins; i++)
        {
    		fprintf(fFIRout,"%g\n", aFilter[i]);
        }
        fclose(fFIRout);
#ifdef ROOT_FOUND
        WriteRootHisto( aFilter, nBins, 0 );
#endif

    }

    void HFSSResponseFileHandlerCore::PrintFIR( fftw_complex* aFilter, int nBins, std::string filename )
    {
    	std::vector<double> vecFilter0;
    	std::vector<double> vecFilter1;
    	LDEBUG( lmclog, "Printing coefficients to file ... ");
    	std::string textFile(filename);
    	std::string suffix(".txt");
    	textFile.replace(textFile.find(".root"),5,suffix);
    	FILE * fFIRout = fopen(textFile.c_str(), "w");
    	for (int i = 0; i < nBins; i++)
    	{
    		fprintf(fFIRout,"%g %g\n", aFilter[i][0], aFilter[i][1]);
    		vecFilter0.push_back(aFilter[i][0]);
    		vecFilter1.push_back(aFilter[i][1]);
    	}
    	fclose(fFIRout);
#ifdef ROOT_FOUND
    	fRootHistoWriter = RootHistoWriter::get_instance();
    	fRootHistoWriter->SetFilename(filename);
    	fRootHistoWriter->OpenFile("RECREATE");
    	fRootHistoWriter->CloseFile();
    	WriteRootHisto( vecFilter0, nBins, 0 );
    	WriteRootHisto( vecFilter1, nBins, 1 );
#endif
    }

    FIRFileHandlerCore::FIRFileHandlerCore():HFSSResponseFileHandlerCore()
    {
    }
    
    FIRFileHandlerCore::~FIRFileHandlerCore()
    {
    }
    
    bool FIRFileHandlerCore::Configure(const scarab::param_node& aParam)
    {
        return true;
    }

    
} /* namespace locust */
