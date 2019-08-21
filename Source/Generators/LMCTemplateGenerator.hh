/*
 * LMC[name]Generator.hh
 *
 *  Created on: Mar 12, 2014
 *      Author: nsoblath
 */

#ifndef LMC[NAME]GENERATOR_HH_
#define LMC[NAME]GENERATOR_HH_

#include "LMCGenerator.hh"

namespace locust
{

    /*!
     @class [name]Generator
     @author N. S. Oblath

     @brief

     @details
     Operates in [domain] space

     Configuration name: "config-name"

     Available configuration options:
     - "param-name": type -- Description

    */
    class [name]Generator : public Generator
    {
        public:
            [name]Generator( const std::string& aName = "config-name" );
            virtual ~[name]Generator();

            bool Configure( const scarab::param_node& aNode );

            void Accept( GeneratorVisitor* aVisitor ) const;

        private:
            bool DoGenerate( Signal* aSignal ) const;

    };

} /* namespace locust */

#endif /* LMC[NAME]GENERATOR_HH_ */
