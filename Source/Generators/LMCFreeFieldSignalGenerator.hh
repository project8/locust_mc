/*
 * LMCFreeFieldSignalGenerator.hh
 *
 *  Created on: Mar 12, 2014
 *      Author: nsoblath
 */

#ifndef LMCFREEFIELDSIGNALGENERATOR_HH_
#define LMCFREEFIELDSIGNALGENERATOR_HH_

#include <KThreeVector.hh>
#include "LMCGenerator.hh"

namespace locust
{

    /*!
     @class FreeFieldSignalGenerator
     @author N. S. Oblath

     @brief

     @details
     Operates in time space

     Configuration name: "kass-signal"

     Available configuration options:
     - "param-name": type -- Description

    */
    class FreeFieldSignalGenerator : public Generator
    {
        public:

            FreeFieldSignalGenerator( const std::string& aName = "freefield-signal" );
            virtual ~FreeFieldSignalGenerator();

            bool Configure( const scarab::param_node* aNode );

            void Accept( GeneratorVisitor* aVisitor ) const;


        private:
            mutable std::vector<KGeoBag::KThreeVector > rReceiver; //Vector that contains 3D position of all points at which the fields are evaluated (ie. along receiver surface)
            mutable std::vector<std::pair<int, double> > PreviousTimes; //Cache the results from previous iteration. [0] is previous index, [1] is corresponding retarded time of previous solution
            double fLO_Frequency;  // typically defined by a parameter in json file.

            std::string gxml_filename;

            mutable std::vector<std::vector<std::array<std::array<double,2>, 3 > > > NFDElectricFieldFreq;  //Should use the KThreeVectors too.....
            mutable std::vector<std::vector<std::array<std::array<double,2>, 3 > > > NFDMagneticFieldFreq;

            bool fWriteNFD;
            std::vector<double> NFDFrequencies;
            std::string fAND_filename;
            std::string fNFD_filename;

            bool DoGenerate( Signal* aSignal ) const;
            void* DriveAntenna(int PreEventCounter, unsigned index, Signal* aSignal) const;

            void NFDWrite() const;

            int FindNode(double tNew, double dtStepSize, int IndexOld) const;

    };

} /* namespace locust */

#endif /* LMCFREEFIELDSIGNALGENERATOR_HH_ */
