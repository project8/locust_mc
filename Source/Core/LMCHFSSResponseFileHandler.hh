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
        virtual bool ReadHFSSFile();
        virtual double ConvolveWithFIRFilter(std::deque<double>);// Convolve input signal (voltage or field) with FIR
        virtual std::pair<double,double> ConvolveWithComplexFIRFilter(std::deque<double> inputBuffer);
        int GetFilterSize() const;//Number of entries in the filter
        double GetFilterResolution() const;//Get the resolution of the filter
        void PrintFIR( std::vector<double>, int nBins, std::string filename );
        void PrintFIR( fftw_complex* aFilter, int nBins, std::string filename );
        bool WriteRootHisto( std::vector<double> aFilter, int nBins, bool bIQ );

        
    protected:
        
        // Member variables
        std::string fHFSSFilename;
        std::vector<double> fFilter;
        fftw_complex* fFilterComplex;
        int fTFNBins;
        int fFIRNBins;
        int fCropIndex;
        double fResolution;
        double fCharacteristicImpedance;
        int fNSkips;
        bool fHFSSFiletype;
        ComplexFFT fComplexFFT;
        bool fIsFIRCreated;
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
    
    inline int HFSSResponseFileHandlerCore::GetFilterSize() const
    {
        return fFIRNBins;
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
        bool ConvertAnalyticTFtoFIR(double initialFreq, std::vector<std::complex<double>> tfArray);
        bool ConvertAnalyticGFtoFIR(std::vector<std::pair<double,std::pair<double,double> > > gfArray);

    
    private:
        //Member variables
        fftw_complex *fTFComplex;
        fftw_complex *fFIRComplex;
        
        // Member functions
        bool ConvertTFtoFIR(std::vector<std::complex<double>> &, bool GeneratedTF);
        bool ConvertStoZ(std::vector<std::complex<double>> &tfArray, bool bConvert);
        bool CropFIR(fftw_complex* anArray, bool bConvert);

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
