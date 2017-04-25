/*
 * LMCKassSignalGenerator.hh
 *
 *  Created on: Mar 12, 2014
 *      Author: nsoblath
 */

#ifndef LMCKASSSIGNALGENERATOR_HH_
#define LMCKASSSIGNALGENERATOR_HH_

#define CENTER_TO_SHORT 0.010 // 10 cm
#define CENTER_TO_ANTENNA 0.010 // 10 cm.
#define PI 3.1415926
#define N_GRID_SIDE 9 //Number of discretized points per side of the receiver
#define LO_FREQUENCY 0.


#include "LMCGenerator.hh"

namespace locust
{

    /*!
     @class KassSignalGenerator
     @author N. S. Oblath

     @brief

     @details
     Operates in time space

     Configuration name: "kass-signal"

     Available configuration options:
     - "param-name": type -- Description

    */
    class KassSignalGenerator : public Generator
    {
        public:

            KassSignalGenerator( const std::string& aName = "kass-signal" );
            virtual ~KassSignalGenerator();

            bool Configure( const ParamNode* aNode );

            void Accept( GeneratorVisitor* aVisitor ) const;


        private:
            std::vector<std::array<double, 3> > rReceiver; //Vector that contains 3D position of all points at which the fields are evaluated (ie. along receiver surface)
            mutable std::vector<std::array<double, 2> > PreviousTimes; //Cache the results from previous iteration. [0] is previous retarded time, [1] is corresponding index
            double fLO_Frequency;  // typically defined by a parameter in json file.

            mutable std::vector<std::vector<std::array<std::array<double,2>, 3 > > > NFDElectricFieldFreq;
            mutable std::vector<std::vector<std::array<std::array<double,2>, 3 > > > NFDMagneticFieldFreq;

            bool fWriteNFD;
            std::vector<double> NFDFrequencies;
            std::string fAND_filename;
            std::string fNFD_filename;

            bool DoGenerate( Signal* aSignal ) const;
            void* DriveAntenna(int PreEventCounter, unsigned index, Signal* aSignal, double *ImaginarySignal) const;
            void* FilterNegativeFrequencies(Signal* aSignal, double *ImaginarySignal) const;

            void NFDWrite() const;

            int FindNode(double tNew, double dtStepSize, int IndexOld) const;
            bool IsInside(double tNew, int IndexMin, int IndexMax) const;
            int BinarySearch(double tNew, int IndexMin, int IndexMax) const;


    };

} /* namespace locust */

#endif /* LMCKASSSIGNALGENERATOR_HH_ */
