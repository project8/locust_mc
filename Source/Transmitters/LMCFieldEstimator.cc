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
	fGeneratorType("FIR"),
	fNFIRFilterBins(-99)
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

    int FieldEstimator::GetFilterSize()
    {
	    return fNFIRFilterBins;
    }

    bool FieldEstimator::ReadFIRFile()
    {
	    fNFIRFilterBins=0;
	    if(!ends_with(fFIRFilename,".txt"))
	    {
		    LERROR(lmclog,"The FIR files should be a .txt file");
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
		    ++fNFIRFilterBins;
	    }
	    fclose(firFile);
	    return true;
    }

    double FieldEstimator::ConvolveWithFIRFilter(Signal *aSignal)
    {	
	double convolution=0.0;
	for(int i=0;i<fNFIRFilterBins;++i)
	{
		// Still needs implementation
		//convolution+=fFIRFilter[i];
	}
	return convolution;
    }

} /* namespace locust */
