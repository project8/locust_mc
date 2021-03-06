/*
 * LMCGeneratorToolbox.hh
 *
 *  Created on: Feb 5, 2014
 *      Author: nsoblath
 */

#ifndef LMCGENERATORTOOLBOX_HH_
#define LMCGENERATORTOOLBOX_HH_
#include "LMCException.hh"

namespace scarab
{
    class param_node;
}

namespace locust
{
    class Generator;

    /*!
     @class GeneratorToolbox
     @author N. S. Oblath

     @brief Creates and configures the requested generators

     @details

     Configuration name:

     Available configuration options:
     - "generators": array -- ordered list of generators

    */
    class GeneratorToolbox
    {
        public:
            GeneratorToolbox();
            virtual ~GeneratorToolbox();

            bool Configure( const scarab::param_node& aNode );

            const Generator* GetFirstGenerator() const;
            Generator* GetFirstGenerator();

        private:
            void ConfigureGenerators( Generator* aGen );

            Generator* fFirstGenerator;

    };

} /* namespace locust */

#endif /* LMCGENERATORTOOLBOX_HH_ */
