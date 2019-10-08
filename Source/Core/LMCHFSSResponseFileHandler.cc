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
    
    // Handler for FIR transmitter
    FIRTransmitterHandler::FIRTransmitterHandler():HFSSResponseFileHandlerCore()
    {
    }
    
    FIRTransmitterHandler::~FIRTransmitterHandler()
    {
    }
    
    bool FIRTransmitterHandler::Configure(const scarab::param_node& aParam)
    {
        if( aParam.has( "fir-transmitter-filename" ) )
        {
            fHFSSFilename=aParam["fir-transmitter-filename"]().as_string();
        }
        if( aParam.has( "fir-transmitter-dt" ) )
        {
            fResolution=aParam["fir-transmitter-dt"]().as_double();
        }
        if( aParam.has( "fir-transmitter-nskips" ) )
        {
            fNSkips=aParam["fir-transmitter-nskips"]().as_int();
        }
        fHFSSFiletype="fir";
        return true;
    }
    
    
    // Handler for FIR receiver
    FIRReceiverHandler::FIRReceiverHandler():HFSSResponseFileHandlerCore()
    {
    }
    
    FIRReceiverHandler::~FIRReceiverHandler()
    {
    }
    
    bool FIRReceiverHandler::Configure(const scarab::param_node& aParam)
    {
        if( aParam.has( "fir-receiver-filename" ) )
        {
            fHFSSFilename=aParam["fir-receiver-filename"]().as_string();
        }
        if( aParam.has( "fir-receiver-dt" ) )
        {
            fResolution=aParam["fir-receiver-dt"]().as_double();
        }
        if( aParam.has( "fir-receiver-nskips" ) )
        {
            fNSkips=aParam["fir-receiver-nskips"]().as_int();
        }
        fHFSSFiletype="fir";
        return true;
    }
    
    // Handler for TF trasnmitter
    TFTransmitterHandler::TFTransmitterHandler():HFSSResponseFileHandlerCore()
    {
    }
    
    TFTransmitterHandler::~TFTransmitterHandler()
    {
    }
    
    bool TFTransmitterHandler::Configure(const scarab::param_node& aParam)
    {
        if( aParam.has( "tf-transmitter-filename" ) )
        {
            fHFSSFilename=aParam["tf-transmitter-filename"]().as_string();
        }
        if( aParam.has( "tf-transmitter-dt" ) )
        {
            fResolution=aParam["tf-transmitter-dt"]().as_double();
        }
        if( aParam.has( "tf-transmitter-nskips" ) )
        {
            fNSkips=aParam["tf-transmitter-nskips"]().as_int();
        }
        fHFSSFiletype="tf";
        return true;
    }
    
    // Handler for TF receiver
    TFReceiverHandler::TFReceiverHandler():HFSSResponseFileHandlerCore()
    {
    }
    
    TFReceiverHandler::~TFReceiverHandler()
    {
    }
    
    bool TFReceiverHandler::Configure(const scarab::param_node& aParam)
    {
        if( aParam.has( "tf-receiver-filename" ) )
        {
            fHFSSFilename=aParam["tf-receiver-filename"]().as_string();
        }
        if( aParam.has( "tf-receiver-dt" ) )
        {
            fResolution=aParam["tf-receiver-dt"]().as_double();
        }
        if( aParam.has( "tf-receiver-nskips" ) )
        {
            fNSkips=aParam["tf-receiver-nskips"]().as_int();
        }
        fHFSSFiletype="tf";
        return true;
    }
    
} /* namespace locust */
