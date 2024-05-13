#ifndef LMCHFSSRESPONSEFILEHANDLER_HH_
#define LMCHFSSRESPONSEFILEHANDLER_HH_

#include <fftw3.h>

#include "param.hh"
#include "LMCComplexFFT.hh"

#ifdef ROOT_FOUND
    #include "LMCRootHistoWriter.hh"
#endif

namespace locust
{
    /*!
     @class HFSSResponseFileHandlerCore
     @author P. T. Surukuchi
     @brief Handles opening of the HFSS-generated files (FIR ot TF), saving them into objects to be used for convolution
     @details
     Available configuration options:
     - "hfss-filetype": string -- The type of file being handler. Currently only Transfer function and Finite Impulse Response
     - "print-fir-debug": bool -- Print text file of FIR and TF coefficients and/or plot Root histograms to file.
     */
    
    class HFSSResponseFileHandlerCore
    {
    public:
        HFSSResponseFileHandlerCore();
        virtual ~HFSSResponseFileHandlerCore();
        
        // Member functions
        virtual bool Configure( const scarab::param_node& aNode);
        virtual bool ReadHFSSFile(int bTE, int l, int m, int n);
        virtual std::pair<double,double> ConvolveWithComplexFIRFilterArray(int bTE, int l, int m, int n, std::deque<double> inputBuffer);
        int GetFilterSizeArray(int bTE, int l, int m, int n) const;//Number of entries in the filter
        int GetNModes() const;
        double GetFilterResolutionArray(int bTE, int l, int m, int n) const;//Get the resolution of the filter
        bool DimensionMultiMode( int nModes );
        void PrintFIR( std::vector<double>, int nBins, std::string filename );
        void PrintFIR( fftw_complex* aFilter, int nBins, std::string filename );
        bool WriteRootHisto( std::vector<double> aFilter, int nBins, bool bIQ );
        double QuadrantCorrection( double aRealValue, double aPhase);

    protected:
        
        // Member variables
        int fNModes;
        int fTFNBins;
        double fResolution;
        std::string fHFSSFilename;
        std::vector< std::vector< std::vector < std::vector< fftw_complex*>>>> fFilterComplexArray;
        std::vector < std::vector < std::vector < std::vector < int >>>> fTFNBinsArray;
        std::vector < std::vector < std::vector < std::vector < int >>>> fFIRNBinsArray;
        std::vector < std::vector < std::vector < std::vector < double >>>> fResolutionArray;
        std::vector < std::vector < std::vector < std::vector < bool >>>> fIsFIRCreatedArray;

        bool fIsFIRCreated;
        double fCharacteristicImpedance;
        int fNSkips;
        bool fHFSSFiletype;
        ComplexFFT fComplexFFT;
        std::string fWindowName;
        double fWindowParam;
        bool fPrintFIR;
        bool fConvertStoZ;
        std::string fOutputPath;

#ifdef ROOT_FOUND
        FileWriter* fRootHistoWriter;
#endif



        //Member functions
        bool ends_with(const std::string &, const std::string &);
    };

    inline int HFSSResponseFileHandlerCore::GetFilterSizeArray(int bTE, int l, int m, int n) const
    {
        return fFIRNBinsArray[bTE][l][m][n];
    }

    inline int HFSSResponseFileHandlerCore::GetNModes() const
    {
        return fNModes;
    }

    inline double HFSSResponseFileHandlerCore::GetFilterResolutionArray(int bTE, int l, int m, int n) const
    {   
        return fResolutionArray[bTE][l][m][n];
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
        bool ReadHFSSFile(int bTE, int l, int m, int n);
        bool ConvertAnalyticTFtoFIR(int bTE, int l, int m, int n, double initialFreq, std::vector<std::complex<double>> tfArray);
        bool ConvertAnalyticGFtoFIR(int nModes, std::vector<std::pair<double,std::pair<double,double> > > gfArray);
    
    private:
        //Member variables
        fftw_complex *fTFComplex;
        fftw_complex *fFIRComplex;
        
        // Member functions
        bool ConvertTFtoFIR(int bTE, int l, int m, int n, std::vector<std::complex<double>> &);
        bool ConvertStoZ(std::vector<std::complex<double>> &tfArray, bool bConvert);

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
    };
} /*namespace locust*/

#endif/*LMCHFSSRESPONSEFILEHANDLER_HH_ */
