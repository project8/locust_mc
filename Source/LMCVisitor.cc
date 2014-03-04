/*
 * LMCVisitor.cc
 *
 *  Created on: Feb 7, 2014
 *      Author: nsoblath
 */

#include "LMCVisitor.hh"

#include "LMCGenerator.hh"
#include "LMCLogger.hh"

namespace locust
{
    LMCLOGGER( lmclog, "GeneratorVisitor" );

    GeneratorVisitor::GeneratorVisitor()
    {
    }

    GeneratorVisitor::~GeneratorVisitor()
    {
    }

    void GeneratorVisitor::Visit( const Generator* aGen )
    {
        LMCWARN( lmclog, "Generic generator-visit function called for generator <" << aGen->GetName() << ">" );
        return;
    }

} /* namespace locust */
