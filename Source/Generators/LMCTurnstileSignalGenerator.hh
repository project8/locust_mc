/*
 * LMCTurnstileSignalGenerator.hh
 *
 *  Created on: Oct 2019
 *      Author: Pranava Teja Surukuchi
 */

#ifndef LMCTURNSTILESIGNALGENERATOR_HH_
#define LMCTURNSTILESIGNALGENERATOR_HH_

#include <fftw3.h>

#include "LMCThreeVector.hh"
#include "LMCGenerator.hh"
#include "LMCChannel.hh"
#include "LMCPatchAntenna.hh"
#include "LMCPowerCombiner.hh"
#include "LMCFieldBuffer.hh"
#include "LMCHilbertTransform.hh"
#include "LMCAntennaSignalTransmitter.hh"
#include "LMCFIRFileHandler.hh"

namespace scarab
{
    class param_node;
}

namespace locust
{
    class Digitizer;
    
    /*!
     @class TurnstileSignalGenerator
     @author  Pranava Teja Surukuchi
     
     @brief Signal generator for a turnstile antenna
     
     @details Uses AntennaSignalTransmitter for obtaining the field which is thewen used to calculate the voltage measured by the patches
     
     Configuration name: "turnstile-signal-generator"
     
     Available configuration options:
     - "input-signal-frequency": double -- Frequency of the incident sine wave.
     - "lo-frequency": double -- Frequency of the local oscillator.
     - "input-signal-amplitude": double -- Amplitude of the incident sine wave.
     - "filter-filename": double -- path to FIR text file.
     - "filter-resolution": double -- time resolution of coefficients in filter-filename.
     - "buffer-size": double -- size of buffer to contain incident E field values.
     - "buffer-margin": double -- distance from beginning of buffer at which to extract Hilbert transform info.  extrapolate to edge across this margin.
     - "domain": string -- Determines whether the sinusoidal test signal is generated in the time
     or frequency domain
     
     Available options: "time"[default] and "freq"
     
     */
    class TurnstileSignalGenerator : public Generator
    {
    public:
        TurnstileSignalGenerator( const std::string& aName = "turnstile-signal-generator" );
        virtual ~TurnstileSignalGenerator();
        
        bool Configure( const scarab::param_node& aNode );
        
        void Accept( GeneratorVisitor* aVisitor ) const;

        double GetRFFrequency() const;
        void SetRFFrequency( double aFrequency );
        
        double GetLOFrequency() const;
        void SetLOFrequency( double aFrequency );
        
        double GetAmplitude() const;
        void SetAmplitude( double aAmplitude );
        
        double GetBufferSize() const;
        void SetBufferSize( double aBufferSize );
        
        double GetBufferMargin() const;
        void SetBufferMargin( double aBufferMargin );
        
        Signal::State GetDomain() const;
        void SetDomain( Signal::State aDomain );
        
        
    private:
        bool DoGenerate( Signal* aSignal );
        bool DoGenerateTime( Signal* aSignal );
        bool DoGenerateFreq( Signal* aSignal );
        bool (TurnstileSignalGenerator::*fDoGenerateFunc)( Signal* aSignal );
        
        double GetAOIFactor(LMCThreeVector TrasmittingPatchNormal,LMCThreeVector ReceivingPatchNormal);
        double RotateZ(int component, double angle, double x, double y);
        bool InitializePatchArray();
        bool InitializePowerCombining();
        double GetVoltageFromField(unsigned channel, unsigned patch,double fieldPhase);
        
        void InitializeBuffers(unsigned filterbuffersize, unsigned fieldbuffersize);
        void FillBuffers(Signal* aSignal, double FieldAmplitude, double FieldPhase, double LOPhase, unsigned index, unsigned channel, unsigned patch);
        void PopBuffers(unsigned channel, unsigned patch);
        void CleanupBuffers();
        
        PowerCombiner fPowerCombiner;
        HilbertTransform fHilbertTransform;
        FIRReceiverHandler fReceiverHandler;
        AntennaSignalTransmitter fAntennaSignalTransmitter;
        std::vector< Channel<PatchAntenna> > allChannels; //Vector that contains pointer to all channels
        std::vector<LMCThreeVector > rReceiver; //Vector that contains 3D position of all points at which the fields are evaluated (ie. along receiver surface)
        double fArrayRadius;  // from json file.
        int fNPatchesPerStrip; // from json file.
        double fPatchSpacing; // from json file.
        std::string gxml_filename;// from json file.
        bool fTextFileWriting;// from json file.
        
        double fRF_frequency;
        double fLO_frequency;
        double fAmplitude;
        unsigned fFieldBufferSize;
        int fSwapFrequency;
        
        std::vector<std::deque<double>> EFieldBuffer;
        std::vector<std::deque<double>> EPhaseBuffer;
        std::vector<std::deque<double>> EAmplitudeBuffer;
        std::vector<std::deque<double>> EFrequencyBuffer;
        std::vector<std::deque<double>> LOPhaseBuffer;
        std::vector<std::deque<unsigned>> IndexBuffer;
        std::vector<std::deque<double>> PatchFIRBuffer;
        fftw_complex *PatchFIRComplex;
    };
    
} /* namespace locust */

#endif /* LMCTURNSTILESIGNALGENERATOR_HH_ */

