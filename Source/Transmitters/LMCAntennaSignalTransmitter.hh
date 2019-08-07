/*
 * LMCAntennaSignalTransmitter.hh
 *
 *  Created on: May 4, 2019
 *      Author: pslocum
 */

#ifndef LMCANTENNASIGNALTRANSMITTER_HH_
#define LMCANTENNASIGNALTRANSMITTER_HH_

#include "LMCThreeVector.hh"
#include "LMCFieldBuffer.hh"
#include "LMCFieldEstimator.hh"
#include "LMCConst.hh"

namespace locust
{

    /*!
     @class AntennaSignalTransmitter
     @author P. L. Slocum

     @brief Class to generate tranmistter from antenna signal

     @details
     Operates in time space

     Configuration name: "antenna-signal"

     Available configuration options:
     - "input-signal-type": 1 -- Leaving an option open for generating different types of signals,
     - "input-signal-frequency": 25.9281e9,
     - "input-signal-amplitude": 1
     - "antenna-x-position": double -- Location of the antenna in the x direction
     - "antenna-y-position": double -- Location of the antenna in the y direction
     - "antenna-z-position": double -- Location of the antenna in the z direction
          
    */
    class AntennaSignalTransmitter
    {
        public:

            AntennaSignalTransmitter();
            virtual ~AntennaSignalTransmitter();

            bool Configure( const scarab::param_node* aNode );

	    /// Generate the electric field based on the voltage input from the config file and convolution with FIR
            double GenerateSignal(Signal *,double acquisitionRate);
	    
	    /// Get initial phase delay 
	    double GetInitialPhaseDelay();
            
	    /// Get the positions of the antenna w.r.t the center of the detector
	    LMCThreeVector GetAntennaPosition() const;
	    
	    /// Set the positions of the antenna w.r.t the center of the detector
	    void SetAntennaPosition(const LMCThreeVector &);
	    
	    /// Initialize the FIR filter and the field estimator 
	    bool InitializeTransmitter();
           
        private:
	    FieldEstimator fFieldEstimator;
	    
	    /// Placeholder for now. Input signal type, 1 for dipole antenna, could be chaged later on.
	    int fInputSignalType;
	    double fInputFrequency;// in GHz
            double fInputAmplitude;// in V/m
	    double fPhaseDelay=0.0; //Delay in the phase that changes for each time sample
	    double fInitialPhaseDelay = 0.0;  //Initial delay in the phase from the the signal arriving from the back of the buffer as well as the delay from signal travel
	    double fArrayRadius=0.0; //Array radius to be used to obtain the phase delay from the tranmistter to the reciever patch
	    
	    double fAntennaPositionX;
	    double fAntennaPositionY;
	    double fAntennaPositionZ;
	    LMCThreeVector fAntennaPosition; // Position of the antenna w.r.t to the center of the array

	    void InitializeBuffers(unsigned); 

	    std::vector<std::deque<double>> delayedVoltageBuffer;

    };

    inline LMCThreeVector AntennaSignalTransmitter::GetAntennaPosition() const
    {
	    return fAntennaPosition;
    } 

    inline void AntennaSignalTransmitter::SetAntennaPosition(const LMCThreeVector &antennaPosition)
    {
	    fAntennaPosition=antennaPosition;
    } 

} /* namespace locust */

#endif /* LMCANTENNASIGNALTRANSMITTER_HH_ */
