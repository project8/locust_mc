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
     - "print-fir-debug": bool -- Print text file of FIR coefficients.
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
        virtual std::pair<double,double> ConvolveWithComplexFIRFilterArray(int l, int m, int n, std::deque<double> inputBuffer);
        virtual std::pair<double,double> ConvolveWithComplexFIRFilter(std::deque<double> inputBuffer);
        int GetFilterSize() const;//Number of entries in the filter
	int GetFilterSizeArray(int l, int m, int n) const;//Number of entries in the filter
        double GetFilterResolution() const;//Get the resolution of the filter
	double GetFilterResolutionArray(int l, int m, int n) const;//Get the resolution of the filter
        void PrintFIR( std::vector<double> );
        void PrintFIR( fftw_complex* aFilter );
        
    protected:
        
        // Member variables
        std::string fHFSSFilename;
        std::vector<double> fFilter;
        fftw_complex* fFilterComplex;
	std::vector< std::vector < std::vector< fftw_complex*>>> fFilterComplexArray;
        int fTFNBins;
        int fFIRNBins;
	std::vector < std::vector < std::vector < int >>> fFIRNBinsArray;
	int fNModes;
        double fResolution;
	std::vector < std::vector < std::vector < double >>> fResolutionArray;
        int fNSkips;
        bool fHFSSFiletype;
        ComplexFFT fComplexFFT;
        bool fIsFIRCreated;
	std::vector < std::vector < std::vector < bool >>> fIsFIRCreatedArray;
        std::string fWindowName;
        double fWindowParam;
        bool fPrintFIR;


        //Member functions
        bool ends_with(const std::string &, const std::string &);
    };
    
    inline int HFSSResponseFileHandlerCore::GetFilterSize() const
    {
        return fFIRNBins;
    }
    
    inline int HFSSResponseFileHandlerCore::GetFilterSizeArray(int l, int m, int n) const
    {
        return fFIRNBinsArray[l][m][n];
    }

    inline double HFSSResponseFileHandlerCore::GetFilterResolution() const
    {
        return fResolution;
    }

    inline double HFSSResponseFileHandlerCore::GetFilterResolutionArray(int l, int m, int n) const
    {   
        return fResolutionArray[l][m][n];
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
        bool ConvertAnalyticTFtoFIR(double initialFreq, std::vector<std::complex<double>> tfArray);
        bool ConvertAnalyticGFtoFIR(int l, int m, int n, std::vector<std::pair<double,std::pair<double,double> > > gfArray);

    
    private:
        //Member variables
        fftw_complex *fTFComplex;
        fftw_complex *fFIRComplex;
        
        // Member functions
        bool ConvertTFtoFIR(std::vector<std::complex<double>> &, bool GeneratedTF);

    protected:
        //Member variables
        double fInitialTFIndex;
        double fTFBinWidth;
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
