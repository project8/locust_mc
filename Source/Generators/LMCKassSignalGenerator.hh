/*
 * LMCKassSignalGenerator.hh
 *
 *  Created on: Mar 12, 2014
 *      Author: nsoblath
 */

#ifndef LMCKASSSIGNALGENERATOR_HH_
#define LMCKASSSIGNALGENERATOR_HH_

#include "LMCGenerator.hh"

#include "LMCKassLocustInterface.hh"
//#include "LMCConst.hh"
//#include "LMCThreeVector.hh"

#include <vector>
using std::vector;


namespace locust
{

    /*!
     @class KassSignalGenerator
     @author N. S. Oblath

     @brief
     Generates a Phase I/ Phase II like signal in waveguide

     @details
     Operates in time space

     Configuration name: "kass-signal"

     Available configuration options:
     - "lo-frequency": double -- Frequency of downmixer
     - "xml-filename": string -- Name of the locust config file

    */
    class KassSignalGenerator : public Generator
    {
        public:

            KassSignalGenerator( const std::string& aName = "kass-signal" );
            virtual ~KassSignalGenerator();

            bool Configure( const scarab::param_node& aNode );

            void Accept( GeneratorVisitor* aVisitor ) const;


        private:
            void KassiopeiaInit(const std::string &aFile);
            void WakeBeforeEvent();
            bool ReceivedKassReady();

            bool DoGenerate( Signal* aSignal );
            void* DriveAntenna(int PreEventCounter, unsigned index, Signal* aSignal, FILE *fp);
            int FindNode(double tNew) const;
            double TE11ModeExcitation() const;
            double TE10ModeExcitation() const;

            double fLO_Frequency;  // typically defined by a parameter in json file.
            bool fTruth; // parameter in json file.  default is false.
            bool fvoltageCheck;  // parameter to print out voltage values.
            std::string gxml_filename;
            std::string gpitchangle_filename;
            double fPhi_t1; // antenna voltage phase in radians.
            double fPhi_t2; // reflecting short voltage phase in radians.
            double fPhiLO_t; // voltage phase of LO in radians;
            int fNPreEventSamples;  // spacing between events.  constant for now, could be randomized.
            mutable double fPreviousRetardedTime;
            mutable int fPreviousRetardedIndex;
            double fEventStartTime;
            bool fEventToFile;
            kl_interface_ptr_t fInterface;


    };

} /* namespace locust */

#endif /* LMCKASSSIGNALGENERATOR_HH_ */
