/*
 * LMCKassSignalGenerator.hh
 *
 *  Created on: Mar 12, 2014
 *      Author: nsoblath
 */

#ifndef LMCKASSSIGNALGENERATOR_HH_
#define LMCKASSSIGNALGENERATOR_HH_

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

            bool Configure( const scarab::param_node* aNode );

            void Accept( GeneratorVisitor* aVisitor ) const;

        private:
            double fLO_Frequency;  // typically defined by a parameter in json file.
            bool DoGenerate( Signal* aSignal ) const;
            void* DriveAntenna(int PreEventCounter, unsigned index, Signal* aSignal) const;
            double TE11ModeExcitation() const;

    };

} /* namespace locust */

#endif /* LMCKASSSIGNALGENERATOR_HH_ */
