/*
 * LMCPatchSignalGenerator.hh
 *
 *  Created on: Feb 19, 2018
 *      Author: pslocum
 */

#ifndef LMCPATCHSIGNALGENERATOR_HH_
#define LMCPATCHSIGNALGENERATOR_HH_

#include "LMCThreeVector.hh"
#include "LMCGenerator.hh"
#include "LMCChannel.hh"
#include "LMCPatchAntenna.hh"
#include "LMCPowerCombiner.hh"
#include "LMCFieldBuffer.hh"
#include "LMCHilbertTransform.hh"


namespace locust
{

    /*!
     @class PatchSignalGenerator
     @author P. L. Slocum

     @brief Generate signal in free space(without wave guide) for phase III and detect with patch array.

     @details
     Operates in time space

     Configuration name: "patch-signal"

     Available configuration options:
     - "param-name": type -- Description
     - "lo-frequency" : double -- the special value tuned down by the local oscillator, e.g., the 24.something giga hertz.
     - "xml-filename" : std::string -- the name of the xml locust config file.
     

    */
    class PatchSignalGenerator : public Generator
    {
        public:

            PatchSignalGenerator( const std::string& aName = "patch-signal" );
            virtual ~PatchSignalGenerator();

            bool Configure( const scarab::param_node& aNode );

            void Accept( GeneratorVisitor* aVisitor ) const;
              
            void AddOnePatchVoltageToStripSum(Signal* aSignal, double VoltageAmplitude, double VoltagePhase, double phi_LO, unsigned channelindex, unsigned z_index, double DopplerFrequency);
            void AddOneFIRVoltageToStripSum(Signal* aSignal, double VoltageFIRSample, double phi_LO, unsigned channelindex, unsigned patchIndex);



        private:
            std::vector< Channel<PatchAntenna> > allChannels; //Vector that contains pointer to all channels
            std::vector<LMCThreeVector > rReceiver; //Vector that contains 3D position of all points at which the fields are evaluated (ie. along receiver surface)
            double fLO_Frequency;  // typically defined by a parameter in json file.
            double fArrayRadius;  // from json file.
            int fNPatchesPerStrip; // from json file.
            double fPatchSpacing; // from json file.
            std::string gxml_filename;
            int fPowerCombiner;
            double fRJunction;
            bool fTextFileWriting;

            double fFilter_resolution;
            std::string gfilter_filename;
            unsigned fFieldBufferSize;
            unsigned fFieldBufferMargin;

            double* GetFIRFilter(int nskips);
            int GetNFilterBins(double* filterarray);
            double GetFIRSample(double* filterarray, int nfilterbins, double dtfilter, unsigned channel, unsigned patch, double AcquisitionRate);
            void InitializeBuffers(unsigned filterbuffersize, unsigned fieldbuffersize);
            void CleanupBuffers();
            void PopBuffers(unsigned channel, unsigned patch);
            void FillBuffers(Signal* aSignal, double DopplerFrequency, double EFieldValue, double LOPhase, unsigned index, unsigned channel, unsigned patch, unsigned dtauConvolutionTime);


            std::vector<std::deque<double>> EFieldBuffer;
            std::vector<std::deque<double>> EPhaseBuffer;
            std::vector<std::deque<double>> EAmplitudeBuffer;
            std::vector<std::deque<double>> EFrequencyBuffer;
            std::vector<std::deque<double>> LOPhaseBuffer;
            std::vector<std::deque<unsigned>> IndexBuffer;
            std::vector<std::deque<double>> PatchFIRBuffer;


            bool DoGenerate( Signal* aSignal );
            void* DriveAntenna(FILE *fp, int PreEventCounter, unsigned index, Signal* aSignal, double* filterarray, unsigned nfilterbins, double dtfilter);
            void InitializePatchArray();

            int FindNode(double tNew) const;
            double GetSpaceTimeInterval(const double &aParticleTime, const double &aReceiverTime, const LMCThreeVector &aParticlePosition, const LMCThreeVector &aReceiverPosition );


            double phiLO_t; // voltage phase of LO in radians;
            double VoltagePhase_t[10000];

    };

} /* namespace locust */

#endif /* LMCPATCHSIGNALGENERATOR_HH_ */
