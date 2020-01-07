#ifndef LMCFIRFILEHANDLER_HH_
#define LMCFIRFILEHANDLER_HH_

#include "LMCException.hh"
#include "LMCHFSSResponseFileHandler.hh"
#include "param.hh"

namespace locust
{
    /*!
     @class FIRTransmitterHandler
     @brief Handles transmitter FIR, saving them into objects to be used for convolution
     @details
     Available configuration options:
     - "fir-transmitter-filename": string -- The location of the file containing impulse response function
     - "fir-transmitter-dt": double (1e-12) -- The size of filter sample width (seconds)
     - "fir-transmitter-nskips": int (1) -- Number of skips to be peformed in reading from the FIR file, this will deteremine the number of FIR bins
     */
    
    class FIRTransmitterHandler: public FIRFileHandlerCore
    {
    public:
        FIRTransmitterHandler();
        virtual ~FIRTransmitterHandler();
        
        // Member functions
        bool Configure( const scarab::param_node& aNode);
    };
    
    /*!
     @class FIRReceiverHandler
     @brief Handles receiver FIR, saving them into objects to be used for convolution
     @details
     Available configuration options:
     - "fir-transmitter-filename": string -- The location of the file containing impulse response function
     - "fir-transmitter-dt": double (1e-12) -- The size of filter sample width (seconds)
     - "fir-transmitter-nskips": int (1) -- Number of skips to be peformed in reading from the FIR file, this will deteremine the number of FIR bins
     */
    
    class FIRReceiverHandler: public FIRFileHandlerCore
    {
    public:
        FIRReceiverHandler();
        virtual ~FIRReceiverHandler();
        
        // Member functions
        bool Configure( const scarab::param_node& aNode);
    };
} /*namespace locust*/

#endif/*LMCFIRFILEHANDLER_HH_ */
