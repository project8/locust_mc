#ifndef LMCHFSSRESPONSEFILEHANDLER_HH_
#define LMCHFSSRESPONSEFILEHANDLER_HH_

#include "LMCException.hh"
#include "param.hh"

namespace locust
{
    /*!
     @class HFSSResponseFileHandlerCore
     @author P. T. Surukuchi
     @brief Handles opening of the HFSS-generated files (FIR ot TF), saving them into objects to be used for convolution
     @details
     Available configuration options:
     - "hfss-filetype": string -- The type of file being handler. Currently only Transfer function and Finite Impulse Response
     */
    
    class HFSSResponseFileHandlerCore
    {
    public:
        HFSSResponseFileHandlerCore();
        virtual ~HFSSResponseFileHandlerCore();
        
        // Member functions
        virtual bool Configure( const scarab::param_node& aNode);
        bool ReadHFSSFile();
        double ConvolveWithFIRFilter(std::deque<double>);// Convolve input signal (voltage or field) with FIR
        int GetFilterSize() const;//Number of entries in the filter
        double GetFilterResolution() const;//Get the resolution of the filter
        
    protected:
        
        // Member variables
        std::string fHFSSFilename;
        std::vector<double> fFilter;
        int fNBins;
        double fResolution;
        int fNSkips;
        bool fHFSSFiletype;
        
        //Member functions
        bool ends_with(const std::string &, const std::string &);
    };
    
    inline int HFSSResponseFileHandlerCore::GetFilterSize() const
    {
        return fNBins;
    }
    
    inline double HFSSResponseFileHandlerCore::GetFilterResolution() const
    {
        return fResolution;
    }
    
    /*!
     @class FIRTransmitterHandler
     @brief Handles transmitter FIR, saving them into objects to be used for convolution
     @details
     Available configuration options:
     - "fir-transmitter-filename": string -- The location of the file containing impulse response function
     - "fir-transmitter-dt": double (1e-12) -- The size of filter sample width (seconds)
     - "fir-transmitter-nskips": int (1) -- Number of skips to be peformed in reading from the FIR file, this will deteremine the number of FIR bins
     */
    
    class FIRTransmitterHandler: public HFSSResponseFileHandlerCore
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
    
    class FIRReceiverHandler: public HFSSResponseFileHandlerCore
    {
    public:
        FIRReceiverHandler();
        virtual ~FIRReceiverHandler();
        
        // Member functions
        bool Configure( const scarab::param_node& aNode);
    };
    
    /*!
     @class TFTransmitterHandler
     @brief Handles transmitter for TF, saving them into objects to be used for convolution
     @details
     Available configuration options:
     - "tf-transmitter-filename": string -- The location of the file containing impulse response function
     - "tf-transmitter-dt": double (1e-12) -- The size of filter sample width (seconds)
     - "tf-transmitter-nskips": int (1) -- Number of skips to be peformed in reading from the FIR file, this will deteremine the number of FIR bins
     */
    
    class TFTransmitterHandler: public HFSSResponseFileHandlerCore
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
    
    class TFReceiverHandler: public HFSSResponseFileHandlerCore
    {
    public:
        TFReceiverHandler();
        virtual ~TFReceiverHandler();
        
        // Member functions
        bool Configure( const scarab::param_node& aNode);
    };
    
} /*namespace locust*/

#endif/*LMCHFSSRESPONSEFILEHANDLER_HH_ */
