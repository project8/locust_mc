/*
 * LMCHFSSFileHandler.cc
 *
 *  Created on: May 11, 2018
 *      Author: P. T. Surukuchi
 */

#include "LMCConst.hh"
#include "LMCHFSSResponseFileHandler.hh"

#include "logger.hh"

namespace locust
{
    LOGGER( lmclog, "HFSSResponseFileHandlerCore" );
    
    HFSSResponseFileHandlerCore::HFSSResponseFileHandlerCore():
    fNBins(1000),
    fResolution(1e-12),
    fNSkips(1),
    fHFSSFiletype("")
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
        fNBins=0;
        if(!ends_with(fHFSSFilename,".txt"))
        {
            LERROR(lmclog,"The FIR file should end in .txt");
            return false;
        }
        double firIndex;
        double filterMagnitude;
        FILE *firFile;
        firFile=fopen(fHFSSFilename.c_str(),"r");
        //logic copied from /LMCPatchSignalGenerator.cc
        int count=0;
        while (!feof(firFile)){
            fscanf(firFile,"%lf %lf",&firIndex,&filterMagnitude);
            if (count%fNSkips==0)
            {
                fFilter.push_back(filterMagnitude);
                ++fNBins;
            }
            ++count;
        }
        fclose(firFile);
        return true;
    }
    
    double HFSSResponseFileHandlerCore::ConvolveWithFIRFilter(std::deque<double> delayedVoltageBuffer)
    {
        double convolution=0.0;
        if(fNBins<=0){
            LERROR(lmclog,"Number of bins in the filter should be positive");
        }
        for(int i=0;i<fNBins;++i)
        {
            convolution+=fFilter[i]*delayedVoltageBuffer[i];
        }
        return convolution;
    }
    
    TFFileHandlerCore::TFFileHandlerCore():HFSSResponseFileHandlerCore()
    {
    }
    
    TFFileHandlerCore::~TFFileHandlerCore()
    {
    }
    
    bool TFFileHandlerCore::Configure(const scarab::param_node& aParam)
    {
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
} /* namespace locust */
