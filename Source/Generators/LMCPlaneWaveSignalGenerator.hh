/*
 * LMCPlaneWaveSignalGenerator.hh
 *
 *  Created on: Feb 9, 2019
 *      Author: pslocum
 */

#ifndef LMCPLANEWAVESIGNALGENERATOR_HH_
#define LMCPLANEWAVESIGNALGENERATOR_HH_

#include "LMCThreeVector.hh"
#include "LMCGenerator.hh"
#include "LMCChannel.hh"
#include "LMCPatchAntenna.hh"
#include "LMCPowerCombiner.hh"
#include "LMCFieldBuffer.hh"
#include "LMCHilbertTransform.hh"
#include "LMCFIRHandler.hh"

namespace locust
{

  /*!
    @class PlaneWaveSignalGenerator
    @author P. L. Slocum

    @brief Generate artibtrary plane wave for calibration and detect with patch array.

    @details
    Operates in time space

    Configuration name: "planewave-signal"

    Available configuration options:
    - "param-name": type -- Description
    - "lo-frequency" : double -- the special value tuned down by the local oscillator, e.g., the 24.something giga hertz.
     
  */
  class PlaneWaveSignalGenerator : public Generator
  {
  public:

    PlaneWaveSignalGenerator( const std::string& aName = "planewave-signal" );
    virtual ~PlaneWaveSignalGenerator();

    bool Configure( const scarab::param_node& aNode );

    void Accept( GeneratorVisitor* aVisitor ) const;
              
//    void AddOnePatchVoltageToStripSum(Signal* aSignal, unsigned bufferIndex, int patchIndex);
    double GetAOIFactor(double AOI, LMCThreeVector PatchNormalVector);
    double GetVoltageAmpFromPlaneWave(int z_index);
    double GetPWPhaseDelayAtPatch(int z_index);

    bool GetPhaseDelay() const;
    void SetPhaseDelay( bool aPhaseDelay );
    bool GetVoltageDamping() const;
    void SetVoltageDamping( bool aVoltageDamping );
    double GetPlaneWaveFrequency() const;
    void SetPlaneWaveFrequency( double aPlaneWaveFrequency );
    double GetLOFrequency() const;
    void SetLOFrequency( double aLOFrequency );
    double GetArrayRadius() const;
    void SetArrayRadius( double aArrayRadius );
    int GetNPatchesPerStrip() const;
    void SetNPatchesPerStrip( int aNPatchesPerStrip );
    double GetPatchSpacing() const;
    void SetPatchSpacing( double aPatchSpacing );
    int GetPowerCombiner() const;
    void SetPowerCombiner( std::string feed );
    double GetAOI() const;
    void SetAOI( double aAOI );
//    bool GetPatchFIRfilter() const;
//    void SetPatchFIRfilter( bool aPatchFIRfilter );
//    std::string GetPatchFIRfilter_filename() const;
//    void SetPatchFIRfilter_filename( std::string aPatchFIRfilterfilename );
//    double GetPatchFIRfilter_resolution() const;
//    void SetPatchFIRfilter_resolution( double aAmplitude );
    double GetAmplitude() const;
    void SetAmplitude( double aAmplitude );
    double GetBufferMargin() const;
    void SetBufferMargin( double aBufferMargin );
    double GetJunctionResistance() const;
    void SetJunctionResistance( double aRJunction );


      
  private:
    // patch and plane wave parameters
    std::vector< Channel<PatchAntenna> > allChannels; //Vector that contains pointer to all channels
    std::vector<LMCThreeVector > rReceiver; //Vector that contains 3D position of all points at which the fields are evaluated (ie. along receiver surface)
    double fLO_Frequency;  // typically defined by a parameter in json file.
    double fRF_Frequency;  // typically defined by a parameter in json file.
    double fArrayRadius;  // from json file.
    int fNPatchesPerStrip; // from json file.
    double fPatchSpacing; // from json file.
    double fAOI; // from json file, in degrees.
    double fRJunction; // for parallel voltage summing

    // options to turn on/off or to select
    int fPowerCombiner; // internally keeps track of power combiner type
    bool fPhaseDelay; // yes/no for calculating phase delays
    bool fVoltageDamping; // yes/no for calculating voltage damping due to junctions
    bool fPatchFIRfilter; // yes/no to use the patch FIR filter
    bool fJunctionCascade; // yes/no to use the S31 junction cascade

//    void ProcessFIRFilter(int nskips);
//    int GetNFilterBins();
    double GetPatchFIRSample(double amp, double startphase, int patchIndex);
//    std::string gpatchfilter_filename;
//    double fPatchFIRfilter_resolution;
//    double FIR_array[1000];
//    int nfilterbins;

      // for FIR filter
    FIRHandler fReceiverFIRHandler;
    double* GetHilbertMagPhase(unsigned bufferIndex);
    double fAmplitude;

    
    bool DoGenerate( Signal* aSignal );
    void* DriveAntenna(int PreEventCounter, unsigned index, Signal* aSignal);
    bool InitializePatchArray();
    bool InitializePowerCombining();
    PowerCombiner testPowerCombiner;

    
    // for buffers
    void InitializeBuffers(unsigned fieldbuffersize);
    void FillBuffers(unsigned bufferIndex, int digitizerIndex, double phiLO, double pwphase, double pwval);
    void PopBuffers(unsigned bufferIndex);

    unsigned fFieldBufferMargin;
    
    std::vector<std::deque<unsigned>> SampleIndexBuffer;
    std::vector<std::deque<double>> LOPhaseBuffer;
    std::vector<std::deque<double>> PWFreqBuffer;
    std::vector<std::deque<double>> PWPhaseBuffer;
    std::vector<std::deque<double>> PWValueBuffer;
    
    std::vector<std::deque<double>> PatchVoltageBuffer;

  };

} /* namespace locust */

#endif /* LMCPLANEWAVESIGNALGENERATOR_HH_ */
