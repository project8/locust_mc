/*
 * LMCRootHistoWriter.cc
 *
 *  Created on: Mar. 21, 2020
 *      Author: pslocum
 */

#include "LMCRootHistoWriter.hh"
#include "logger.hh"



namespace locust
{
    LOGGER( lmclog, "RootHistoWriter" );


    RootHistoWriter::RootHistoWriter()
    {
    }

    RootHistoWriter::~RootHistoWriter()
    {
    }

    bool RootHistoWriter::Configure( const scarab::param_node& aParam )
    {


    	if( !FileWriter::Configure(aParam))
    	{
    		LERROR(lmclog,"Error configuring FileWriter class from RootHistoWriter child class");
    	}

    	return true;
    }

    void RootHistoWriter::Write1DHisto(TH1D* aHisto)
    {
    	aHisto->Write();
    }

    void RootHistoWriter::Write2DHisto(TH2D* aHisto)
    {
    	aHisto->Write();
    }

    void RootHistoWriter::WriteVector1DHisto(std::vector<double> aVector, double xmin, double xmax)
    {
    	TH1D* aHisto = new TH1D("histo", "title", aVector.size(), xmin, xmax);
    	unsigned count = 1;
    	for (auto it = aVector.begin(); it!=aVector.end(); ++it)
    	{
    	    aHisto->SetBinContent(count,*it);
    	    count ++;
    	}
    	aHisto->Write();
    }

} /* namespace locust */
