/*
 * LMCWGNoiseSourceGenerator.hh
 *
 *  Created on: Mar 12, 2014
 *      Author: nsoblath
 */

#ifndef LMCWGNOISESOURCEGENERATOR_HH_
#define LMCWGNOISESOURCEGENERATOR_HH_

#include "LMCGenerator.hh"

#include <complex>
#include <vector>

namespace locust
{

    /*!
     @class WGNoiseSourceGenerator
     @author N. S. Oblath

     @brief Simulates the waveguide noise sources

     @details
     Operates in frequency space

     Configuration name: "noise-source"

     Available configuration options:
     - "param-name": type -- Description

    */
    class WGNoiseSourceGenerator : public Generator
    {
        public:
            WGNoiseSourceGenerator( const std::string& aName = "noise-source" );
            virtual ~WGNoiseSourceGenerator();

            bool Configure( const scarab::param_node* aNode );

            void Accept( GeneratorVisitor* aVisitor ) const;

        private:
            bool DoGenerate( Signal* aSignal ) const;

            void GenerateWGe( unsigned aFreqSize, double aFreqBW );

            std::vector< std::complex< double > > fWGe;

    };

} /* namespace locust */

#endif /* LMCWGNOISESOURCEGENERATOR_HH_ */
