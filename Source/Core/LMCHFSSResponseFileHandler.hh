#ifndef LMCHFSSRESPONSEFILEHANDLER_HH_
#define LMCHFSSRESPONSEFILEHANDLER_HH_

#include <fftw3.h>

#include "param.hh"
#include "LMCComplexFFT.hh"

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
        virtual bool ReadHFSSFile();
        virtual double ConvolveWithFIRFilter(std::deque<double>);// Convolve input signal (voltage or field) with FIR
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
        virtual bool Configure( const scarab::param_node& aNode) override;
        bool ReadHFSSFile() override;
    
    private:
        //Member variables
        fftw_complex *fTFComplex;
        fftw_complex *fFIRComplex;
        ComplexFFT fIFFT;
        
        //Member functions
        bool ConvertTFtoFIR(std::vector<std::complex<double>> &);
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
        virtual bool Configure( const scarab::param_node& aNode) override;
        bool ReadHFSSFile() override;
    };
} /*namespace locust*/

#endif/*LMCHFSSRESPONSEFILEHANDLER_HH_ */
