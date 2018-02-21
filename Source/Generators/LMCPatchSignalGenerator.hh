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

            bool Configure( const scarab::param_node* aNode );

            void Accept( GeneratorVisitor* aVisitor ) const;


        private:
            std::vector<LMCThreeVector > rReceiver; //Vector that contains 3D position of all points at which the fields are evaluated (ie. along receiver surface)
            std::vector<std::pair<int, double> > PreviousTimes; //Cache the results from previous iteration. [0] is previous index, [1] is corresponding retarded time of previous solution
            double fLO_Frequency;  // typically defined by a parameter in json file.

            std::string gxml_filename;

            std::vector<std::vector<std::array<std::array<double,2>, 3 > > > NFDElectricFieldFreq;  //Should use the LMCThreeVectors too.....
            std::vector<std::vector<std::array<std::array<double,2>, 3 > > > NFDMagneticFieldFreq;

            bool fWriteNFD;
            std::vector<double> NFDFrequencies;
            std::string fAND_filename;
            std::string fNFD_filename;

            bool DoGenerate( Signal* aSignal );
            void* DriveAntenna(int PreEventCounter, unsigned index, Signal* aSignal);

            void NFDWrite();

            int FindNode(double tNew, double dtStepSize, int IndexOld) const;
            double m(double angle, double x0, double y0);
            double b(double m, double x0, double y0);
            double directivity(double m1, double m2, double x0, double y0, double x, double y);
            double GetSpaceTimeInterval(const double &aParticleTime, const double &aReceiverTime, const LMCThreeVector &aParticlePosition, const LMCThreeVector &aReceiverPosition );





    };

} /* namespace locust */

#endif /* LMCPATCHSIGNALGENERATOR_HH_ */
