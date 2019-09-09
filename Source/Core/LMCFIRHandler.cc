/*
 * LMCFIRHandler.cc
 *
 *  Created on: May 11, 2018
 *      Author: P. T. Surukuchi
 */

#include "LMCConst.hh"
#include "LMCFIRHandler.hh"

#include "logger.hh"

namespace locust
{
    LOGGER( lmclog, "FIRHandler" );
    
    FIRHandler::FIRHandler():
    fFIRFilename("blank.txt"),
    fNFIRFilterBins(-99),
    fFilterResolution(1e-12)
    {
    }
    
    FIRHandler::~FIRHandler()
    {
    }
    
    bool FIRHandler::Configure(const scarab::param_node& aParam)
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
    
    bool FIRHandler::ends_with(const std::string &str, const std::string &suffix)
    {
        //copied from https://stackoverflow.com/a/20446239
        return str.size() >= suffix.size() &&
        str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
    }
    
    bool FIRHandler::ReadFIRFile()
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
    
    int FIRHandler::GetFilterSize()
    {
        return fNFIRFilterBins;
    }
    
    double FIRHandler::GetFilterResolution()
    {
        return fFilterResolution;
    }
    
    double FIRHandler::ConvolveWithFIRFilter(std::deque<double> delayedVoltageBuffer)
    {
        double convolution=0.0;
        for(int i=0;i<fNFIRFilterBins;++i)
        {
            convolution+=fFIRFilter[i]*delayedVoltageBuffer[i];
        }
        return convolution;
    }
    
} /* namespace locust */
