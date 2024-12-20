/*
 * LMCArraySignalGenerator.hh
 *
 *  Created on: Feb 19, 2018
 *      Author: pslocum
 */

#ifndef LMCARRAYSIGNALGENERATOR_HH_
#define LMCARRAYSIGNALGENERATOR_HH_

#include "LMCPowerCombinerInterface.hh"
#include "LMCTransmitterInterface.hh"
#include "LMCAntennaElementPositionerInterface.hh"
#include "LMCThreeVector.hh"
#include "LMCGenerator.hh"
#include "LMCChannel.hh"
#include "LMCFieldBuffer.hh"
#include "LMCHilbertTransform.hh"
#include "LMCFIRFileHandler.hh"
#include "LMCTFFileHandler.hh"
#include "LMCKassLocustInterface.hh"
#include <vector>
#include <algorithm>    // std::min
#include "LMCException.hh"


namespace locust
{

    /*!
     @class ArraySignalGenerator
     @author P. L. Slocum

     @brief Generate signal in free space(without wave guide) for phase III and detect with patch array.

     @details
     Operates in time space

     Configuration name: "array-signal"

     Available configuration options:
     - "param-name": type -- Description
     - "lo-frequency" : double -- local oscillator frequency
     - "xml-filename" : std::string -- the name of the xml locust config file.
     - "buffer-size" :  std::int -- number of elements in deque buffers to contain field information.
     	 	 These buffers control arrival times of fields and provide a short time series for the Hilbert transform.
     - "lo-frequency":  local oscillator frequency in Hz.
     - "array-radius":  radius of cylindrical antenna array in meters.
     - "npatches-per-strip":  number of patch antennas on each strip.
     - "patch-spacing":  spacing between patches on one strip in meters.
     - "zshift-array":  shift of whole antenna array along z axis, for testing (meters).
     - "swap-frequency":  number of digitizer samples after which buffer memory is reset.  This
     	 	 becomes more important for large numbers of patches
     - "allow-fast-sampling": use sampling interval to define overlap time between incident field and FIR.
     	 	 default is false, which sets overlap to be the entire duration of FIR.

    */

    class ArraySignalGenerator : public Generator
    {
        public:

            ArraySignalGenerator( const std::string& aName = "array-signal" );
            virtual ~ArraySignalGenerator();

            bool Configure( const scarab::param_node& aNode );
            bool ConfigureInterface(Signal* aSignal);
            bool RecordRunParameters(Signal* aSignal);
            const scarab::param_node* GetParameters();
            void SetParameters( const scarab::param_node& aNode );


            void Accept( GeneratorVisitor* aVisitor ) const;
              



        private:
            std::vector< Channel<Receiver*> > allRxChannels; //Vector of channels with pointers to Rx elements.
            double fLO_Frequency;
            int fNPreEventSamples;  // spacing between events.  constant for now, could be randomized.
            int fThreadCheckTime;  // time (ms) to check for response from Kass thread.
            double fArrayRadius;
            int fNElementsPerStrip;
            int fNSubarrays;
            double fZShiftArray;
            double fElementSpacing;
            std::string gxml_filename;
            bool fTextFileWriting;
            unsigned fFieldBufferSize;
            unsigned fFIRzeroBuffer;
            int fSwapFrequency;
            bool fKassNeverStarted;
            bool fSkippedSamples;
            bool fAllowFastSampling;
            double fphiLO; // voltage phase of LO in radians;

            void KassiopeiaInit(const std::string &aFile);
            void WakeBeforeEvent();
            bool TryWakeAgain();
            bool ReceivedKassReady();

        	void InitializeFieldPoints(std::vector< Channel<Receiver*> > allRxChannels);
            double GetFIRSample(int nFilterBinsRequired, double dtfilter, unsigned channel, unsigned patch);
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
            std::vector<std::deque<double>> ElementFIRBuffer;
            std::vector<std::deque<double>> FIRfrequencyBuffer;


            bool DoGenerate( Signal* aSignal );
            bool DriveAntenna(int startingIndex, unsigned index, Signal* aSignal, int nFilterBinsRequired, double dtfilter);
            bool InitializeElementArray();
            AntennaElementPositioner* fAntennaElementPositioner;
            Transmitter* fTransmitter; // transmitter object
            PowerCombiner* fPowerCombiner;
            TFReceiverHandler fTFReceiverHandler;
            HilbertTransform fHilbertTransform;

            kl_interface_ptr_t fInterface;
            const scarab::param_node* fParam;



    };

} /* namespace locust */

#endif /* LMCARRAYSIGNALGENERATOR_HH_ */
