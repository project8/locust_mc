/*
 * LMCCavitySignalGenerator.hh
 *
 *  Created on: Mar 30, 2021
 *      Author: pslocum
 */

#ifndef LMCCAVITYSIGNALGENERATOR_HH_
#define LMCCAVITYSIGNALGENERATOR_HH_

#include "LMCThreeVector.hh"
#include <boost/math/special_functions/bessel.hpp>
#include "LMCGenerator.hh"
#include "LMCChannel.hh"
#include "LMCKassTransmitter.hh"
#include "LMCKassLocustInterface.hh"
#include "LMCFieldBuffer.hh"
#include <vector>
#include <sstream>
#include <string>
#include "LMCException.hh"
#include <algorithm>    // std::min



namespace locust
{

    /*!
     @class CavitySignalGenerator
     @author P. L. Slocum

     @brief Generate signal in cavity and detect it with antennas.

     @details
     Operates in [tbd] domain

     Configuration name: "cavity-signal"

     Available configuration options:
     - "param-name": type -- Description
     - "lo-frequency" : double -- local oscillator frequency
     - "xml-filename" : std::string -- the name of the xml locust config file.
     - "lo-frequency":  local oscillator frequency in Hz.
     - "array-radius":  radius of cylindrical antenna array in meters.
     - "nelements-per-strip":  number of patch antennas on each strip.
     - "element-spacing":  spacing between patches on one strip in meters.
     - "zshift-array":  shift of whole antenna array along z axis, for testing (meters).

    */

    class CavitySignalGenerator : public Generator
    {
        public:

            CavitySignalGenerator( const std::string& aName = "cavity-signal" );
            virtual ~CavitySignalGenerator();

            bool Configure( const scarab::param_node& aNode );

            void Accept( GeneratorVisitor* aVisitor ) const;
              
            Signal::State GetDomain() const;
            void SetDomain( Signal::State aDomain );
            std::vector<double> TE_E(int l, int m, int n, double r, double theta, double z) const;
            std::vector<double> TE_H(int l, int m, int n, double r, double theta, double z) const;
            std::vector<double> TM_E(int l, int m, int n, double r, double theta, double z) const;
            std::vector<double> TM_H(int l, int m, int n, double r, double theta, double z) const;
            double Integrate(int l, int m, int n, bool teMode, bool eField);
            void CheckNormalization();
            void PrintModeMaps();




        private:
//            std::vector< Channel<Receiver*> > allRxChannels; //Vector of channels with pointers to Rx elements.
            double fR, fL; // cavity radius, length
            int fnPixels;
            double fLO_Frequency;
            int fNModes;
            int fNPreEventSamples;  // spacing between events.  constant for now, could be randomized.
            int fThreadCheckTime;  // time (ms) to check for response from Kass thread.
            double fArrayRadius;
            std::string gxml_filename;
            bool fKassNeverStarted;
            bool fSkippedSamples;
            double fphiLO; // voltage phase of LO in radians;



            void ReadBesselZeroes(std::string filename, std::vector<std::vector<double> > &data);
            void KassiopeiaInit(const std::string &aFile);
            void WakeBeforeEvent();
            bool ReceivedKassReady();

//        	void InitializeFieldPoints(std::vector< Channel<Receiver*> > allRxChannels);

            bool DoGenerate( Signal* aSignal );
            bool DoGenerateTime( Signal* aSignal );
            bool DoGenerateTimeKass( Signal* aSignal );
            bool DoGenerateFreq( Signal* aSignal );
            bool (CavitySignalGenerator::*fDoGenerateFunc)( Signal* aSignal );

            Transmitter* fTransmitter; // transmitter object

            kl_interface_ptr_t fInterface;


    };

} /* namespace locust */

#endif /* LMCCAVITYSIGNALGENERATOR_HH_ */
