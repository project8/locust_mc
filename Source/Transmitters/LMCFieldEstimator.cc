/*
 * LMCFieldEstimator.cc
 *
 *  Created on: May 11, 2018
 *      Author: P. T. Surukuchi
 */

#include "LMCConst.hh"
#include "LMCFieldEstimator.hh"

#include "logger.hh"

namespace locust
{
    LOGGER( lmclog, "FieldEstimator" );

    FieldEstimator::FieldEstimator():
	fFIRFilename("blank.txt"),
	fNFIRFilterBins(-99),
	fFilterResolution(1e-12)
    {
    }

    FieldEstimator::~FieldEstimator()
    {
    }

    bool FieldEstimator::Configure(const scarab::param_node& aParam)
    {
	    if( aParam.has( "fir-filename" ) )
	    {
		    fFIRFilename=aParam["fir-filename"]().as_string();
	    }
	    if( aParam.has( "filter-dt" ) )
	    {
		    fFilterResolution=aParam["filter-dt"]().as_double();
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

    int FieldEstimator::GetFilterSize()
    {
	    return fNFIRFilterBins;
    }

    double FieldEstimator::GetFilterResolution()
    {
	    return fFilterResolution;
    }

    double FieldEstimator::ConvolveWithFIRFilter(std::deque<double> delayedVoltageBuffer)
    {	
	double convolution=0.0;
	for(int i=0;i<fNFIRFilterBins;++i)
	{
		convolution+=fFIRFilter[i]*delayedVoltageBuffer[i];
	}
	return convolution;
    }

    double FieldEstimator::ApplyDerivative(double voltagePhase)
    {
	    return -sin(voltagePhase);
    }

    double FieldEstimator::GetFieldAtOrigin(double inputAmplitude,double voltagePhase)
    {
	    //double normalizedVoltage = cos(voltagePhase);
	    double normalizedDerivative = ApplyDerivative(voltagePhase);
	    // Only missing tau, f_g
	    // And distance will be applied later
	    double field = inputAmplitude*normalizedDerivative/(2*LMCConst::Pi()*LMCConst::C());
	    return field; 
    }
    

} /* namespace locust */
