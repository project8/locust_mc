#ifndef LMCFIRHANDLER_HH_
#define LMCFIRHANDLER_HH_

#include "LMCException.hh"
#include "param.hh"

namespace locust
{
    /*!
     @class FIRHandler
     @author P. T. Surukuchi
     @brief Handles opening of the FIR files, saving them into objects to be used for convolution
     @details
     Available configuration options:
     - "fir-filename": string -- The location of the file containing impulse response
     - "filter-dt": double (1e-12) -- The size of filter sample width (seconds)
     */
    
    class FIRHandler
    {
    public:
        FIRHandler();
        virtual ~FIRHandler();
        
        // Member functions
        bool Configure( const scarab::param_node& aNode );
        bool ReadFIRFile();
        double ConvolveWithFIRFilter(std::deque<double>);// Convolve input signal (voltage or field) with FIR
        int GetFilterSize();//Number of entries in the filter
        double GetFilterResolution();//Get the resolution of the filter
        
    private:
        
        // Member variables
        std::string fFIRFilename;
        std::vector<double> fFIRFilter;
        int fNFIRFilterBins;
        double fFilterResolution;
        
        //Member functions
        bool ends_with(const std::string &, const std::string &);
    };
    
} /*namespace locust*/

#endif/*LMCFIRHANDLER_HH_ */
