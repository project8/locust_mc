/*
 * LMCFIRFileHandler.cc
 *
 *  Created on: May 11, 2018
 *      Author: P. T. Surukuchi
 */

#include "LMCConst.hh"
#include "LMCFIRFileHandler.hh"
#include "logger.hh"

namespace locust
{
    LOGGER( lmclog, "FIRFileHandlerCore" );
    
    // Handler for FIR transmitter
    FIRTransmitterHandler::FIRTransmitterHandler():FIRFileHandlerCore()
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
    FIRReceiverHandler::FIRReceiverHandler():FIRFileHandlerCore()
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
    
} /* namespace locust */
