/*
 * LMCHFSSFileHandler.cc
 *
 *  Created on: May 11, 2018
 *      Author: P. T. Surukuchi
 */

#include "LMCConst.hh"
#include "LMCException.hh"
#include "LMCHFSSResponseFileHandler.hh"

#include "logger.hh"

#include<unistd.h>

namespace locust
{
    LOGGER( lmclog, "HFSSResponseFileHandlerCore" );
    
    HFSSResponseFileHandlerCore::HFSSResponseFileHandlerCore():
    fTFNBins(1000),
    fFIRNBins(2000),
    fResolution(1e-12),
    fNSkips(1),
    fComplexFFT(),
    fHFSSFiletype(""),
    fIsFIRCreated(false)
    {
    }
    
    HFSSResponseFileHandlerCore::~HFSSResponseFileHandlerCore()
    {
    }
    
    bool HFSSResponseFileHandlerCore::Configure(const scarab::param_node& aParam)
    {
        return true;
    }
    
    bool HFSSResponseFileHandlerCore::ends_with(const std::string &str, const std::string &suffix)
    {
        //copied from https://stackoverflow.com/a/20446239
        return str.size() >= suffix.size() &&
        str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
    }
    
    bool HFSSResponseFileHandlerCore::ReadHFSSFile()
    {
        return true;
    }
    
    double HFSSResponseFileHandlerCore::ConvolveWithFIRFilter(std::deque<double> inputBuffer)
    {
        double convolution=0.0;
        if(fFIRNBins<=0){
            LERROR(lmclog,"Number of bins in the filter should be positive");
        }
	int firBinNumber=0;
	for (auto it = inputBuffer.begin();it!=inputBuffer.end(); ++it)
	{
	    convolution+=*(it)*fFilter[firBinNumber];
	    firBinNumber++;
	}

        return convolution;
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
    }
    
    bool TFFileHandlerCore::Configure(const scarab::param_node& aParam)
    {
        return true;
    }
  
    bool TFFileHandlerCore::ConvertTFtoFIR(std::vector<std::complex<double>> &tfArray)
    {
        if(fTFNBins<=0)
        {
            LERROR(lmclog,"The size of transfer function has to be positive integer");
            return false;
        }
        //Might need to be moved to a different function
        fTFComplex=(fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fTFNBins);
        
        for (int i = 0; i < fTFNBins; ++i)
        {
            fTFComplex[i][0]=tfArray.at(i).real();
            fTFComplex[i][1]=tfArray.at(i).imag();
        }
        fComplexFFT.SetupIFFT(fTFNBins,fInitialTFIndex,fTFBinWidth);
	fFIRNBins=fTFNBins+2*fComplexFFT.GetShiftNBins();
        fFIRComplex=(fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fFIRNBins);
        fComplexFFT.GenerateFIR(fTFNBins,fTFComplex,fFIRComplex);
	fResolution=fComplexFFT.GetTimeResolution();
        for (int i = 0; i < fFIRNBins; ++i){
            fFilter.push_back(fFIRComplex[i][0]);
        }
        LDEBUG( lmclog, "Finished IFFT to convert transfer function to FIR");
        return true;
    }
    
    bool TFFileHandlerCore::ReadHFSSFile()
    {
	if(fIsFIRCreated) 
	{
	    return true;	
	}
        fTFNBins=0;
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
        //logic copied from /LMCPatchSignalGenerator.cc
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
		    if(fTFNBins==0)fInitialTFIndex=tfIndex*pow(10.0,9);
                    const std::complex<double> temp(tfRealValue,tfImaginaryValue);
                    tfArray.push_back(temp);
                    ++fTFNBins;
                }
            }
        }
        tfFile.close();
        LDEBUG( lmclog, "Finished reading transfer function file");
        if(!ConvertTFtoFIR(tfArray)){
            return false;
        }
	fIsFIRCreated=true;
        return true;
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
    
    bool FIRFileHandlerCore::ReadHFSSFile()
    {
	if(fIsFIRCreated) 
	{
	    return true;	
	}
        fFIRNBins=0;
        if(!ends_with(fHFSSFilename,".txt"))
        {
            LERROR(lmclog,"The FIR file " << fHFSSFilename.c_str() <<" doesn't end in .txt");
            return false;
        }
        double firIndex;
        double filterMagnitude;
        FILE *firFile;
        firFile=fopen(fHFSSFilename.c_str(),"r");

	if(!access(fHFSSFilename.c_str(), F_OK ))
	{
            LERROR(lmclog,"The FIR file " << fHFSSFilename.c_str() <<" doesn't exist");
            return false;
	}
        //logic copied from /LMCPatchSignalGenerator.cc
        int count=0;
        
        while (!feof(firFile)){
            fscanf(firFile,"%lf %lf",&firIndex,&filterMagnitude);
            if (count%fNSkips==0)
            {
                fFilter.push_back(filterMagnitude);
                ++fFIRNBins;
            }
            ++count;
        }
        fclose(firFile);
        LDEBUG( lmclog, "Finished reading FIR file");
	fIsFIRCreated=true;
        return true;
    }
    
} /* namespace locust */
