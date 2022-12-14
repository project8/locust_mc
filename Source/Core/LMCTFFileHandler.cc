/*
 * LMCHFSSFileHandler.cc
 *
 *  Created on: May 11, 2018
 *      Author: P. T. Surukuchi
 */

#include "LMCConst.hh"
#include "LMCTFFileHandler.hh"

#include "logger.hh"

namespace locust
{
    LOGGER( lmclog, "TFFileHandlerCore" );
    // Handler for TF trasnmitter
    TFTransmitterHandler::TFTransmitterHandler():TFFileHandlerCore()
    {
    }
    
    TFTransmitterHandler::~TFTransmitterHandler()
    {
    }
    
    bool TFTransmitterHandler::Configure(const scarab::param_node& aParam)
    {
        if(!fComplexFFT.Configure(aParam))
        {
            LERROR(lmclog,"Error configuring ComplexFFT class");
        }
        if( aParam.has( "tf-transmitter-filename" ) )
        {
            fHFSSFilename=aParam["tf-transmitter-filename"]().as_string();
        }
        if( aParam.has( "tf-transmitter-bin-width" ) )
        {
            fTFBinWidth=aParam["tf-transmitter-bin-width"]().as_double();
        }
        if( aParam.has( "tf-transmitter-nskips" ) )
        {
            fNSkips=aParam["tf-transmitter-nskips"]().as_int();
        }
        if(aParam.has("tf-window-type"))
        {
            fWindowName = aParam["tf-window-type"]().as_string();
        }
        if(aParam.has("tf-window-parameter"))
        {
            fWindowParam = aParam["tf-window-parameter"]().as_double();
        }
        fHFSSFiletype="tf";
        return true;
    }
    
    // Handler for TF receiver
    TFReceiverHandler::TFReceiverHandler():TFFileHandlerCore()
    {
    }
    
    TFReceiverHandler::~TFReceiverHandler()
    {
    }
    
    bool TFReceiverHandler::Configure(const scarab::param_node& aParam)
    {


    	if( !HFSSResponseFileHandlerCore::Configure(aParam))
    	{
    		LERROR(lmclog,"Error configuring HFSSResponseFileHandler class from TFReceiverHandler subclass");
    	}

        if(!fComplexFFT.Configure(aParam))
        {
            LERROR(lmclog,"Error configuring ComplexFFT class");
        }
        if( aParam.has( "tf-receiver-filename" ) )
        {
            fHFSSFilename=aParam["tf-receiver-filename"]().as_string();
        }
        
        if( aParam.has( "tf-receiver-bin-width" ) )
        {
            fTFBinWidth=aParam["tf-receiver-bin-width"]().as_double();
        }
        if( aParam.has( "tf-receiver-nskips" ) )
        {
            fNSkips=aParam["tf-receiver-nskips"]().as_int();
        }
        if(aParam.has("tf-window-type"))
        {
            fWindowName = aParam["tf-window-type"]().as_string();
        }
        if(aParam.has("tf-window-parameter"))
        {
            fWindowParam = aParam["tf-window-parameter"]().as_double();
        }
        fHFSSFiletype="tf";
        return true;
    }
    
} /* namespace locust */
