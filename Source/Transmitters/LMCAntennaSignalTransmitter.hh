/*
 * LMCAntennaSignalTransmitter.hh
 *
 *  Created on: May 4, 2019
 *      Author: pslocum
 */

#ifndef LMCANTENNASIGNALTRANSMITTER_HH_
#define LMCANTENNASIGNALTRANSMITTER_HH_

#include "LMCTransmitter.hh"
#include "LMCTransmitterHardware.hh"
#include "LMCDipoleAntenna.hh"
#include "LMCTurnstileAntenna.hh"
#include "LMCThreeVector.hh"
#include "LMCFieldBuffer.hh"
#include "LMCFIRFileHandler.hh"
#include "LMCTFFileHandler.hh"
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
     - "transmitter-frequency": 0.0,
     - "antenna-voltage-amplitude": 1
     - "transmitter-antenna-type":  string "antenna-signal-dipole" or "antenna-signal-turnstile"
     
     */
    class AntennaSignalTransmitter : public Transmitter
    {
    public:
        
        AntennaSignalTransmitter();
        virtual ~AntennaSignalTransmitter();
        
        bool Configure( const scarab::param_node& aNode );
        
        /// Generate the electric field based on the voltage input from the config file and convolution with FIR
        virtual double* GetEFieldCoPol(int fieldPointIndex, double dt);
        
        /// Get initial phase delay
        double GetInitialPhaseDelay();
        
        /// Initialize the FIR filter and the field estimator
        bool InitializeTransmitter();
        
        /// Select dipole or turnstile
        bool SetAntennaType( std::string antennaType );
	    
    private:
        TFTransmitterHandler fTransmitterHandler;
        TransmitterHardware* fTransmitterHardware;
        
        /// Placeholder for now. Input signal type, 1 for dipole antenna, could be chaged later on.
        int fInputSignalType;
        double fInputFrequency;// in GHz
        double fInputAmplitude;// in V/m
        double fPhaseDelay=0.0; //Delay in the phase that changes for each time sample
        double fInitialPhaseDelay = 0.0;  //Initial delay in the phase from the the signal arriving from the back of the buffer as well as the delay from signal travel
        int fAntennaType;
	    
	//Add incidentKVector  
        void AddIncidentKVector(LMCThreeVector pointOfInterest);
        
        void AddPropagationPhaseDelay(LMCThreeVector pointOfInterest);

	//Apply derivative of a given signal. This will be more complicated with implmentation of other field types
        //PTS: Move this to a core file sometime later
        double ApplyDerivative(double voltagePhase);
        
        //Get the value of the field at the transmitter origin for a given amplitude and phase.
        double GetFieldAtOrigin(double inputAmplitude,double voltagePhase);

        void InitializeBuffers(unsigned);
        
        std::vector<std::deque<double>> delayedVoltageBuffer;
        
    };
    
    inline double AntennaSignalTransmitter::ApplyDerivative(double voltagePhase)
    {
        return -sin(voltagePhase);
    }
    
} /* namespace locust */

#endif /* LMCANTENNASIGNALTRANSMITTER_HH_ */
