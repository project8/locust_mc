/*
 * LMCFieldParameterGenerator.cc
 *
 *  Created on: Mar 04, 2020
 *      Author: P. T. Surukuchi
 */

#include "LMCFieldParameterGenerator.hh"

namespace locust
{
    LOGGER( lmclog, "FieldParameterGenerator" );

    MT_REGISTER_GENERATOR(FieldParameterGenerator, "field-parameter");

    FieldParameterGenerator::FieldParameterGenerator( const std::string& aName ):
        FreeSpaceGenerator( aName )
    {
    }

    FieldParameterGenerator::~FieldParameterGenerator()
    {
    }

    bool FieldParameterGenerator::Configure( const scarab::param_node& aParam )
    {
	    FreeSpaceGenerator::Configure(aParam);
    }
} /* namespace locust */

