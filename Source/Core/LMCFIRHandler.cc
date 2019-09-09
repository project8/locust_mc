/*
 * LMCFIRHandler.cc
 *
 *  Created on: May 11, 2018
 *      Author: P. T. Surukuchi
 */

#include "LMCConst.hh"
#include "LMCFIRHandler.hh"

#include "logger.hh"

namespace locust
{
    LOGGER( lmclog, "FIRHandler" );
    
    FIRHandler::FIRHandler():
    fNFIRBins(-99),
    fFIRResolution(1e-12),
    fNFIRSkips(1)
    {
    }
    
    FIRHandler::~FIRHandler()
    {
    }
    
    bool FIRHandler::Configure(const scarab::param_node& aParam)
    {
        if( aParam.has( "fir-filename" ) )
        {
            fFIRFilename=aParam["fir-filename"]().as_string();
        }
        if( aParam.has( "fir-dt" ) )
        {
            fFIRResolution=aParam["fir-dt"]().as_double();
        }
        if( aParam.has( "fir-nskips" ) )
        {
            fNFIRSkips=aParam["fir-nskips"]().as_int();
        }
        return true;
    }
    
    bool FIRHandler::ends_with(const std::string &str, const std::string &suffix)
    {
        //copied from https://stackoverflow.com/a/20446239
        return str.size() >= suffix.size() &&
        str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
    }
    
    bool FIRHandler::ReadFIRFile()
    {
        fNFIRBins=0;
        if(!ends_with(fFIRFilename,".txt"))
        {
            LERROR(lmclog,"The FIR file should end in .txt");
            return false;
        }
        double firIndex;
        double filterMagnitude;
        FILE *firFile;
        firFile=fopen(fFIRFilename.c_str(),"r");
        //logic copied from /LMCPatchSignalGenerator.cc
        int count=0;
        while (!feof(firFile)){
            fscanf(firFile,"%lf %lf",&firIndex,&filterMagnitude);
            if (count%fNFIRSkips==0)
            {
                fFIRFilter.push_back(filterMagnitude);
                ++fNFIRBins;
            }
            ++count;
        }
        fclose(firFile);
        return true;
    }
    
    double FIRHandler::ConvolveWithFIRFilter(std::deque<double> delayedVoltageBuffer)
    {
        double convolution=0.0;
        if(fNFIRBins<=0){
            LERROR(lmclog,"Number of bins in the filter should be positive");
        }
        for(int i=0;i<fNFIRBins;++i)
        {
            convolution+=fFIRFilter[i]*delayedVoltageBuffer[i];
        }
        return convolution;
    }
    
} /* namespace locust */
