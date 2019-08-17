#ifndef LMCFIELDESTIMATOR_HH_
#define LMCFIELDESTIMATOR_HH_ 

//#include "LMCThreeVector.hh"
#include "LMCGenerator.hh"
#include "LMCException.hh"
#include "param.hh"

/*
#include <string>
#include <vector>
#include <array>
#include <sstream>
*/
namespace locust
{
    /*!
         @class FieldEstimator
         @author P. T. Surukuchi
         @brief Estimates field based on a specified voltage provided to a particular antenna
         @details Estimates field based on a specified voltage provided to a particular antenna. Currently does this for a given HFSS-generated impulse response. Includes the convolution of FIR if it is used for the estimation of fields. 
	 Available configuration options:
     	 - "generator-type": string -- Define if the generator to be used is based on impulse response or a analytical
	 - "fir-filename": string -- The location of the file containing impulse response
	 - "filter-dt": double (1e-12) -- The size of filter sample width (seconds)
    */
    
    class FieldEstimator 
    {
        public:
            FieldEstimator();
            virtual ~FieldEstimator();
	    
	    // Member functions
	    bool Configure( const scarab::param_node& aNode );
	    bool ReadFIRFile();
	    double ConvolveWithFIRFilter(std::deque<double>);// Convolve input signal with FIR 
	    int GetFilterSize();//Number of entries in the filter
	    double GetFilterResolution();//Get the resolution of the filter
	    //Apply derivative of a given signal. This will be more complicated with implmentation of other field types 
	    double ApplyDerivative(double voltagePhase);
	    //Get the value of the field for a given amplitude and phase. 
	    //Perhaps has to be moved to LMCAntennaSignalTransmitter
	    double GetFieldAtOrigin(double inputAmplitude,double voltagePhase); 

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

#endif/*LMCFIELDESTIMATOR_HH_ */
