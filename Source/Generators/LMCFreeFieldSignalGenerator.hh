/*
 * LMCFreeFieldSignalGenerator.hh
 *
 *  Created on: Mar 12, 2014
 *      Author: nsoblath
 */

#ifndef LMCFREEFIELDSIGNALGENERATOR_HH_
#define LMCFREEFIELDSIGNALGENERATOR_HH_

#include "LMCThreeVector.hh"
#include "LMCGenerator.hh"

#include "LMCChannel.hh"
#include "LMCPatchAntenna.hh"
#include "LMCLienardWiechert.hh"
#include "LMCKassLocustInterface.hh"


namespace locust
{

    /*!
     @class FreeFieldSignalGenerator
     @author N. S. Oblath

     @brief Generate signal in free space(without wave guide) for phase III

     @details
     Operates in time space

     Configuration name: "freefield-signal"

     Available configuration options:
     - "param-name": type -- Description
     - "lo-frequency" : double -- the special value tuned down by the local oscillator, e.g., the 24.something giga hertz.
     - "xml-filename" : std::string -- the name of the xml locust config file.
     - "and-filename" : std::string -- the file of the hfss config file.
     

    */
    class FreeFieldSignalGenerator : public Generator
    {
        public:

            FreeFieldSignalGenerator( const std::string& aName = "freefield-signal" );
            virtual ~FreeFieldSignalGenerator();

            bool Configure( const scarab::param_node& aNode );

            void Accept( GeneratorVisitor* aVisitor ) const;


        private:
            std::vector< Channel<PatchAntenna> > allChannels; //Vector that contains pointer to all channels

            double fLO_Frequency;  // typically defined by a parameter in json file.
            int fNPreEventSamples;  // spacing between events.  constant for now, could be randomized.
            double fArrayRadius;
            int fNElementsPerStrip;
            double fElementSpacing;
            bool fCorporateFeed;
            bool fPileupMode; //simulate tracks in sequence or as piling up?
            int fPileupSeed; 
            LienardWiechert fFieldSolver;

            kl_interface_ptr_t fInterface;

            std::string gxml_filename;

            bool DoGenerate( Signal* aSignal );
            void* DriveAntenna(int PreEventCounter, unsigned index, Signal* aSignal);
            void InitializePatchArray();
            void KassiopeiaInit(const std::string &aFile);
            bool WakeBeforeEvent();
            bool ReceivedKassReady();

    };

} /* namespace locust */

#endif /* LMCFREEFIELDSIGNALGENERATOR_HH_ */
