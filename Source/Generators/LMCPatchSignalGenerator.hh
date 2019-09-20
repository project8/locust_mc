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
#include "LMCFIRHandler.hh"


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
     - "lo-frequency" : double -- local oscillator frequency
     - "xml-filename" : std::string -- the name of the xml locust config file.
     

    */
    class PatchSignalGenerator : public Generator
    {
        public:

            PatchSignalGenerator( const std::string& aName = "patch-signal" );
            virtual ~PatchSignalGenerator();

            bool Configure( const scarab::param_node& aNode );

            void Accept( GeneratorVisitor* aVisitor ) const;
              



        private:
            std::vector< Channel<PatchAntenna> > allChannels; //Vector that contains pointer to all channels
            double fLO_Frequency;  // typically defined by a parameter in json file.
            double fArrayRadius;  // from json file.
            int fNPatchesPerStrip; // from json file.
            double fPatchSpacing; // from json file.
            std::string gxml_filename;
            bool fTextFileWriting;
            unsigned fFieldBufferSize;
            double fphiLO; // voltage phase of LO in radians;

            bool WakeBeforeEvent();
            bool ReceivedKassReady();
            double GetAOIFactor(LMCThreeVector IncidentKVector, double PatchPhi);
            double GetEFieldCoPol(PatchAntenna* currentPatch, LMCThreeVector IncidentElectricField, LMCThreeVector IncidentKVector, double PatchPhi, double DopplerFrequency);
            void RecordIncidentFields(FILE *fp, LMCThreeVector IncidentMagneticField, LMCThreeVector IncidentElectricField, LMCThreeVector IncidentKVector, double PatchPhi, double DopplerFrequency);
            double GetFIRSample(int nfilterbins, double dtfilter, unsigned channel, unsigned patch);
            void InitializeBuffers(unsigned filterbuffersize, unsigned fieldbuffersize);
            void CleanupBuffers();
            void PopBuffers(unsigned channel, unsigned patch);
            void FillBuffers(Signal* aSignal, double DopplerFrequency, double EFieldValue, double LOPhase, unsigned index, unsigned channel, unsigned patch);


            std::vector<std::deque<double>> EFieldBuffer;
            std::vector<std::deque<double>> EPhaseBuffer;
            std::vector<std::deque<double>> EAmplitudeBuffer;
            std::vector<std::deque<double>> EFrequencyBuffer;
            std::vector<std::deque<double>> LOPhaseBuffer;
            std::vector<std::deque<unsigned>> IndexBuffer;
            std::vector<std::deque<double>> PatchFIRBuffer;


            bool DoGenerate( Signal* aSignal );
            void DriveAntenna(FILE *fp, int PreEventCounter, unsigned index, Signal* aSignal, int nfilterbins, double dtfilter);
            bool InitializePatchArray();
            bool InitializePowerCombining();
            FIRHandler fReceiverFIRHandler;
            PowerCombiner fPowerCombiner;
            HilbertTransform fHilbertTransform;

            int FindNode(double tNew) const;
            double GetSpaceTimeInterval(const double &aParticleTime, const double &aReceiverTime, const LMCThreeVector &aParticlePosition, const LMCThreeVector &aReceiverPosition );

    };

} /* namespace locust */

#endif /* LMCPATCHSIGNALGENERATOR_HH_ */
