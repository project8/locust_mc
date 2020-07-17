/*
 * LMCRootGraphWriter.cc
 *
 *  Created on: Mar. 21, 2020
 *      Author: pslocum
 */

#include "LMCRootGraphWriter.hh"
#include "logger.hh"



namespace locust
{
    LOGGER( lmclog, "RootGraphWriter" );


    RootGraphWriter::RootGraphWriter()
    {
    }

    RootGraphWriter::~RootGraphWriter()
    {
    }

    bool RootGraphWriter::Configure( const scarab::param_node& aParam )
    {

    	if( !FileWriter::Configure(aParam))
    	{
    		LERROR(lmclog,"Error configuring FileWriter class from RootGraphWriter child class");
    	}

    	return true;
    }

    void RootGraphWriter::Write2DGraph(TGraph* aGraph)
    {
    	aGraph->Write();
    }

    void RootGraphWriter::WriteVectorGraph(std::vector<double> xVector, std::vector<double> yVector)
    {
    	int n = xVector.size();
    	double* xArray = new double[n];
    	double* yArray = new double [n];

    	unsigned count = 0;
    	for (auto it = xVector.begin(); it!=xVector.end(); ++it)
    	{
    	    xArray[count] = *it;
    	    count ++;
    	}
    	count = 0;
    	for (auto it = yVector.begin(); it!=yVector.end(); ++it)
    	{
    	    yArray[count] = *it;
    	    count ++;
    	}

    	TGraph* aGraph = new TGraph(n, xArray, yArray);
    	aGraph->Write();


    }



} /* namespace locust */
