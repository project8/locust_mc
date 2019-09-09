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
     - "fir-transmitter(receiver)-filename": string -- The location of the file containing impulse response
     - "fir-transmitter(receiver)-dt": double (1e-12) -- The size of filter sample width (seconds)
     - "fir-transmitter(receiver)-nskips": int (1) -- Number of skips to be peformed in reading from the FIR file, this will deteremine the number of FIR bins
     */
    
    class FIRHandler
    {
    public:
        FIRHandler();
        virtual ~FIRHandler();
        
        // Member functions
        bool Configure( const scarab::param_node& aNode, bool isTransmitter=false);
        bool ReadFIRFile();
        double ConvolveWithFIRFilter(std::deque<double>);// Convolve input signal (voltage or field) with FIR
        int GetFilterSize() const;//Number of entries in the filter
        double GetFilterResolution() const;//Get the resolution of the filter
        
    private:
        
        // Member variables
        std::string fFIRFilename;
        std::vector<double> fFIRFilter;
        int fNFIRBins;
        double fFIRResolution;
        int fNFIRSkips;
        bool fIsTransmitter;//Ad-hoc definition of FIR filter type if true transmitter other-wise receiver
        
        //Member functions
        //Check weather the given input string ends with another string
        //PTS: Should be moved into a core class
        bool ends_with(const std::string &, const std::string &);
    };
    
    inline int FIRHandler::GetFilterSize() const
    {
        return fNFIRBins;
    }
    
    inline double FIRHandler::GetFilterResolution() const
    {
        return fFIRResolution;
    }
} /*namespace locust*/

#endif/*LMCFIRHANDLER_HH_ */
