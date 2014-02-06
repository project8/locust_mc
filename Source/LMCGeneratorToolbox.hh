/*
 * LMCGeneratorToolbox.hh
 *
 *  Created on: Feb 5, 2014
 *      Author: nsoblath
 */

#ifndef LMCGENERATORTOOLBOX_HH_
#define LMCGENERATORTOOLBOX_HH_

#include "RandomLib/Random.hpp"

namespace locust
{
    class Generator;
    class ParamNode;

    class GeneratorToolbox
    {
        public:
            GeneratorToolbox();
            virtual ~GeneratorToolbox();

            void Configure( const ParamNode* aNode );

            const Generator* GetFirstGenerator() const;

        private:
            void ConfigureGenerators( Generator* aGen );

            Generator* fFirstGenerator;

            RandomLib::Random fRNG;

    };

} /* namespace locust */

#endif /* LMCGENERATORTOOLBOX_HH_ */
