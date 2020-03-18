/*
 * LMCFieldParameterGenerator.hh
 *
 *  Created on: March 04, 2020
 *      Author: P. T. Surukuchi
 */

#ifndef LMCFIELDPARAMETERGENERATOR_HH_
#define LMCFIELDPARAMETERGENERATOR_HH_

#include "LMCTransmitterInterfaceGenerator.hh"

namespace locust
{

    /*!
     @class FieldParameterGenerator
     @author P. T. Surukuchi

     @brief Generate signal in free space(without wave guide) for phase III and estimate field parameters (Poynting vector etc) 

     @details
     Operates in time space

     Configuration name: "field-parameter"

     Available configuration options:

    */

    class FieldParameterGenerator : public TransmitterInterfaceGenerator 
    {
        public:

            FieldParameterGenerator( const std::string& aName = "field-parameter" );
            virtual ~FieldParameterGenerator();

            bool Configure( const scarab::param_node& aNode );
    };

} /* namespace locust */

#endif /* LMCFIELDPARAMETERGENERATOR_HH_ */
