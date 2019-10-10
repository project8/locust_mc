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
     @class TFFileHandlerCore
     @brief Handles TF file handler for core class, saving them into objects to be used for convolution
     @details
     Available configuration options:
     */
    
    class TFFileHandlerCore: public HFSSResponseFileHandlerCore
    {
    public:
        TFFileHandlerCore();
        virtual ~TFFileHandlerCore();
        
        // Member functions
        virtual bool Configure( const scarab::param_node& aNode);
    };
    
    /*!
     @class TFHandlerCore
     @brief Handles TF file handler for core class, saving them into objects to be used for convolution
     @details
     Available configuration options:
     */
    
    class FIRFileHandlerCore: public HFSSResponseFileHandlerCore
    {
    public:
        FIRFileHandlerCore();
        virtual ~FIRFileHandlerCore();
        
        // Member functions
        virtual bool Configure( const scarab::param_node& aNode);
    };
} /*namespace locust*/

#endif/*LMCHFSSRESPONSEFILEHANDLER_HH_ */
