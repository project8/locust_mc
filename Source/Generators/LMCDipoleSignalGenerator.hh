/*
 * LMCDipoleSignalGenerator.hh
 *
 *  Created on: July 11 2019
 *      Author: Pranava Teja Surukuchi
 */

#ifndef LMCDIPOLESIGNALGENERATOR_HH_
#define LMCDIPOLESIGNALGENERATOR_HH_ 

#include "LMCThreeVector.hh"
#include "LMCGenerator.hh"
#include "LMCChannel.hh"
#include "LMCPatchAntenna.hh"
#include "LMCPowerCombiner.hh"
#include "LMCFieldBuffer.hh"
#include "LMCHilbertTransform.hh"
#include "LMCAntennaSignalTransmitter.hh"


namespace scarab
{
  class param_node;
}

namespace locust
{
  class Digitizer;

    /*!
     @class DipoleSignalGenerator
     @author  Pranava Teja Surukuchi 

     @brief 

     @details

     Configuration name: "dipole-signal-generator"

     Available configuration options:
     - "input-signal-frequency": double -- Frequency of the incident sine wave.
     - "lo-frequency": double -- Frequency of the local oscillator.
     - "amplitude": double -- Amplitude of the incident sine wave.
     - "filter-filename": double -- path to FIR text file.
     - "filter-resolution": double -- time resolution of coefficients in filter-filename.
     - "buffer-size": double -- size of buffer to contain incident E field values.
     - "buffer-margin": double -- distance from beginning of buffer at which to extract Hilbert transform info.  extrapolate to edge across this margin.
     - "domain": string -- Determines whether the sinusoidal test signal is generated in the time 
            or frequency domain
    
     Available options: "time" and "freq" [default]

    */
    class DipoleSignalGenerator : public Generator
    {
        public:
            DipoleSignalGenerator( const std::string& aName = "antenna-dipole-generator" );
            virtual ~DipoleSignalGenerator();

            bool Configure( const scarab::param_node* aNode );

            void Accept( GeneratorVisitor* aVisitor ) const;

	    void AddOnePatchVoltageToStripSum(Signal* aSignal, double VoltageAmplitude, double VoltagePhase, double phi_LO, unsigned channelindex, unsigned z_index, double DopplerFrequency);

	    void AddOneFIRVoltageToStripSum(Signal* aSignal, double VoltageFIRSample, double phi_LO, unsigned channelindex, unsigned patchIndex);

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
            bool (DipoleSignalGenerator::*fDoGenerateFunc)( Signal* aSignal );

	    void* DriveAntenna(FILE *fp, int PreEventCounter, unsigned index, Signal* aSignal, double* filterarray, unsigned nfilterbins, double dtfilter);
	    double RotateZ(int component, double angle, double x, double y);
	    void InitializePatchArray();

            double* GetFIRFilter(int nskips);
            int GetNFilterBins(double* filterarray);
            double GetFIRSample(double* filterarray, int nfilterbins, double dtfilter, unsigned channel, unsigned patch, double AcquisitionRate);

            void InitializeBuffers(unsigned filterbuffersize, unsigned fieldbuffersize);
            void FillBuffers(Signal* aSignal, double FieldAmplitude, double FieldPhase, double LOPhase, unsigned index, unsigned channel, unsigned patch, unsigned dtauConvolutionTime);
            void PopBuffers(unsigned channel, unsigned patch);
            void CleanupBuffers();

            AntennaSignalTransmitter fAntennaSignalGenerator; 
	    std::vector< Channel<PatchAntenna> > allChannels; //Vector that contains pointer to all channels
	    std::vector<LMCThreeVector > rReceiver; //Vector that contains 3D position of all points at which the fields are evaluated (ie. along receiver surface)
            double fArrayRadius;  // from json file.
            int fNPatchesPerStrip; // from json file.
            double fPatchSpacing; // from json file.
            std::string gxml_filename;
            int fPowerCombiner;
            bool fTextFileWriting;

	    double* filterarray;
            double fRF_frequency;
            double fLO_frequency;
            double fAmplitude;
            double fFilter_resolution;
            std::string gfilter_filename;
            unsigned fFieldBufferSize;
            unsigned fFieldBufferMargin;
            unsigned fNPatches; 

            std::vector<std::deque<double>> EFieldBuffer;
            std::vector<std::deque<double>> EPhaseBuffer;
            std::vector<std::deque<double>> EAmplitudeBuffer;
            std::vector<std::deque<double>> EFrequencyBuffer;
            std::vector<std::deque<double>> LOPhaseBuffer;
            std::vector<std::deque<unsigned>> IndexBuffer;
            std::vector<std::deque<double>> PatchFIRBuffer;
            std::vector<std::deque<unsigned>> ConvolutionTimeBuffer;
    };

} /* namespace locust */

#endif /* LMCDIPOLESIGNALGENERATOR_HH_ */
