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

    bool Configure( const scarab::param_node* aNode );

    void Accept( GeneratorVisitor* aVisitor ) const;
              
    void AddOnePatchVoltageToStripSum(Signal* aSignal, double VoltageAmplitude, double VoltagePhase, double phi_LO, unsigned channelindex, unsigned z_index, double DopplerFrequency);
    double GetAOIFactor(double AOI, LMCThreeVector PatchNormalVector);
    double GetVoltageAmpFromPlaneWave(int z_index);
    double GetTimeDelayToPatch(int z_index);
      
  private:
    std::vector< Channel<PatchAntenna> > allChannels; //Vector that contains pointer to all channels
    std::vector<LMCThreeVector > rReceiver; //Vector that contains 3D position of all points at which the fields are evaluated (ie. along receiver surface)
    double fLO_Frequency;  // typically defined by a parameter in json file.
    double fRF_Frequency;  // typically defined by a parameter in json file.
    double fArrayRadius;  // from json file.
    int fNPatchesPerStrip; // from json file.
    double fPatchSpacing; // from json file.
    int fPowerCombiner; // internally keeps track of power combiner type
    bool fPhaseDelay; // yes/no for calculating phase delays
    bool fVoltageDamping; // yes/no for calculating voltage damping due to junctions
    double fAOI; // from json file, in degrees.

    // for FIR filter 
    void ProcessFIRFilter(int nskips);
    int GetNFilterBins();
    double GetPatchFIRSample(int bufferIndex);
    bool fPatchFIRfilter;
    std::string gpatchfilter_filename;
    double fPatchFIRfilter_resolution;
    double fAmplitude;
    double FIR_array[1000];
    int nfilterbins;
    
    bool DoGenerate( Signal* aSignal );
    void* DriveAntenna(int PreEventCounter, unsigned index, Signal* aSignal);
    void InitializePatchArray();

    double phiLO_t; // voltage phase of LO in radians;

    void InitializeBuffers(unsigned fieldbuffersize);
    void FillPWBuffers(double timeSampleSize, int digitizerIndex, int channelIndex, int patchIndex, int bufferIndex);
    void PopBuffers(unsigned bufferIndex);
    
    std::vector<std::vector<double>> SampleIndexBuffer;
    std::vector<std::vector<double>> LOPhaseBuffer;
    std::vector<std::vector<double>> PWPhaseBuffer;
    std::vector<std::vector<double>> PWMagBuffer;
    std::vector<std::vector<double>> PatchVoltageBuffer;

  };

} /* namespace locust */

#endif /* LMCPLANEWAVESIGNALGENERATOR_HH_ */
