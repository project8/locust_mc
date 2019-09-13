/*
 * LMCFIRHandlerCore.cc
 *
 *  Created on: May 11, 2018
 *      Author: P. T. Surukuchi
 */

#include "LMCConst.hh"
#include "LMCFIRHandler.hh"

#include "logger.hh"

namespace locust
{
    LOGGER( lmclog, "FIRHandlerCore" );
    
    FIRHandlerCore::FIRHandlerCore():
    fNFIRBins(-99),
    fFIRResolution(1e-12),
    fNFIRSkips(1)
    {
    }
    
    FIRHandlerCore::~FIRHandlerCore()
    {
    }
    
    bool FIRHandlerCore::Configure(const scarab::param_node& aParam)
    {
        return true;
    }
    
    bool FIRHandlerCore::ends_with(const std::string &str, const std::string &suffix)
    {
        //copied from https://stackoverflow.com/a/20446239
        return str.size() >= suffix.size() &&
        str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
    }
    
    bool FIRHandlerCore::ReadFIRFile()
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
    
    double FIRHandlerCore::ConvolveWithFIRFilter(std::deque<double> delayedVoltageBuffer)
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
    
    FIRTransmitterHandler::FIRTransmitterHandler():FIRHandlerCore()
    {
    }
    
    FIRTransmitterHandler::~FIRTransmitterHandler()
    {
    }
    
    bool FIRTransmitterHandler::Configure(const scarab::param_node& aParam)
    {
        if( aParam.has( "fir-transmitter-filename" ) )
        {
            fFIRFilename=aParam["fir-transmitter-filename"]().as_string();
        }
        if( aParam.has( "fir-transmitter-dt" ) )
        {
            fFIRResolution=aParam["fir-transmitter-dt"]().as_double();
        }
        if( aParam.has( "fir-transmitter-nskips" ) )
        {
            fNFIRSkips=aParam["fir-transmitter-nskips"]().as_int();
        }
        return true;
    }
    
    FIRReceiverHandler::FIRReceiverHandler():FIRHandlerCore()
    {
    }
    
    FIRReceiverHandler::~FIRReceiverHandler()
    {
    }
    
    bool FIRReceiverHandler::Configure(const scarab::param_node& aParam)
    {
        if( aParam.has( "fir-receiver-filename" ) )
        {
            fFIRFilename=aParam["fir-receiver-filename"]().as_string();
        }
        if( aParam.has( "fir-receiver-dt" ) )
        {
            fFIRResolution=aParam["fir-receiver-dt"]().as_double();
        }
        if( aParam.has( "fir-receiver-nskips" ) )
        {
            fNFIRSkips=aParam["fir-receiver-nskips"]().as_int();
        }
        return true;
    }
    
} /* namespace locust */
