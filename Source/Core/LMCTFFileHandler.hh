#ifndef LMCTFFILEHANDLER_HH_
#define LMCTFFILEHANDLER_HH_

#include "LMCException.hh"
#include "LMCHFSSResponseFileHandler.hh"
#include "param.hh"

namespace locust
{
    /*!
     @class TFTransmitterHandler
     @brief Handles transmitter for TF, saving them into objects to be used for convolution
     @details
     Available configuration options:
     - "tf-transmitter-filename": string -- The location of the file containing impulse response function
     - "tf-transmitter-dt": double (1e-12) -- The size of filter sample width (seconds)
     - "tf-transmitter-nskips": int (1) -- Number of skips to be peformed in reading from the FIR file, this will deteremine the number of FIR bins
     */
    
    class TFTransmitterHandler: public TFFileHandlerCore
    {
    public:
        TFTransmitterHandler();
        virtual ~TFTransmitterHandler();
        
        // Member functions
        bool Configure( const scarab::param_node& aNode);
    };
    
    
    /*!
     @class TFReceiverHandler
     @brief Handles receiver for transfer functions, saving them into objects to be used for convolution
     @details
     Available configuration options:
     - "tf-transmitter-filename": string -- The location of the file containing impulse response function
     - "tf-transmitter-dt": double (1e-12) -- The size of filter sample width (seconds)
     - "tf-transmitter-nskips": int (1) -- Number of skips to be peformed in reading from the FIR file, this will deteremine the number of FIR bins
     */
    
    class TFReceiverHandler: public TFFileHandlerCore
    {
    public:
        TFReceiverHandler();
        virtual ~TFReceiverHandler();
        
        // Member functions
        bool Configure( const scarab::param_node& aNode);
    };
    
} /*namespace locust*/

#endif/*LMCTFFILEHANDLER_HH_ */
