/*
 * LMCFieldEstimator.cc
 *
 *  Created on: May 11, 2018
 *      Author: P. T. Surukuchi
 */

#include "LMCConst.hh"
#include "LMCFieldEstimator.hh"

/*#include <algorithm>
#include <iostream>
#include <fstream>
#include <regex>
#include <math.h>

#include <stdlib.h>
#include <time.h>
*/
#include "logger.hh"

namespace locust
{
    LOGGER( lmclog, "FieldEstimator" );

    FieldEstimator::FieldEstimator():
	fGeneratorType("FIR")
    {
    }

    FieldEstimator::~FieldEstimator()
    {
    }

    bool FieldEstimator::Configure(const scarab::param_node* aParam)
    {
	    if( aParam == NULL) return true;
	    if( aParam->has( "fir-filename" ) )
	    {
		    fFIRFilename=aParam->get_value<std::string>("fir-filename");
	    }
	    return true;
    }

    bool FieldEstimator::ends_with(const std::string &str, const std::string &suffix)
    { 
	    //copied from https://stackoverflow.com/a/20446239
	    return str.size() >= suffix.size() &&
			           str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
    }

    bool FieldEstimator::ReadFIRFile()
    {
	    if(!ends_with(fFIRFilename,".txt"))
	    {
		    LERROR(lmclog,"The FIR files should be a .csv file");
		    return false;
	    }
	    double firIndex;
	    double filterMagnitude;
	    FILE *firFile;
	    firFile=fopen(fFIRFilename.c_str(),"r");
	    //logic copied from /LMCPatchSignalGenerator.cc
	    while (!feof(firFile)){
		    fscanf(firFile,"%lf %lf",&firIndex,&filterMagnitude);
		    fFIRFilter.push_back(filterMagnitude);
	    }
	    fclose(firFile);
	    return true;
    }

} /* namespace locust */
