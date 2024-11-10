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
        int fNModes = 2;
        if( aParam.has( "n-modes" ) )
        {   
            fNModes = aParam["n-modes"]().as_int();
        }  
        if( aParam.has( "fir-transmitter-filename" ) )
        {
            fHFSSFilename=aParam["fir-transmitter-filename"]().as_string();
        }
        if( aParam.has( "fir-transmitter-dt" ) )
        {
	    for(int bTE = 0; bTE<2; bTE++)
	    {
		for(int l = 0; l<fNModes; l++)
		{
			for(int m = 0; m<fNModes; m++)
			{
				for(int n = 0; n<fNModes; n++)
				{
            				fResolutionArray[bTE][l][m][n]=aParam["fir-transmitter-dt"]().as_double(); //Should eventually be implemented to read in specific values for different modes
				}
			}
		}
	    }
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
        int fNModes = 2;
        if( aParam.has( "n-modes" ) ) 
        {   
                fNModes = aParam["n-modes"]().as_int();
        }   
        if( aParam.has( "fir-receiver-filename" ) )
        {
            fHFSSFilename=aParam["fir-receiver-filename"]().as_string();
        }
        if( aParam.has( "fir-receiver-dt" ) )
        {
            for(int bTE = 0; bTE<2; bTE++)
            {   
                for(int l = 0; l<fNModes; l++)
                {   
                        for(int m = 0; m<fNModes; m++)
                        {   
                                for(int n = 0; n<fNModes; n++)
                                {   
                                        fResolutionArray[bTE][l][m][n]=aParam["fir-receiver-dt"]().as_double(); //Should eventually be implemented to read in specific values for different modes
                                }   
                        }   
                }   
            }   
        }
        if( aParam.has( "fir-receiver-nskips" ) )
        {
            fNSkips=aParam["fir-receiver-nskips"]().as_int();
        }
        fHFSSFiletype="fir";
        return true;
    }
    
} /* namespace locust */
