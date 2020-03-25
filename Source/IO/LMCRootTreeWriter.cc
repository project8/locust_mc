/*
 * LMCRootTreeWriter.cc
 *
 *  Created on: Mar. 21, 2020
 *      Author: pslocum
 */

#include "LMCRootTreeWriter.hh"
#include "logger.hh"
using std::string;



namespace locust
{
    LOGGER( lmclog, "RootTreeWriter" );


    RootTreeWriter::RootTreeWriter():
	fFile ( 0 )
    {
    }

    RootTreeWriter::~RootTreeWriter()
    {
    }

    bool RootTreeWriter::Configure( const scarab::param_node& aParam )
    {

/*    	if( !FileWriter::Configure(aParam))
    	{
    		LERROR(lmclog,"Error configuring FileWriter class from RootTreeWriter child class");
    	}
*/
    	return true;
    }

    double RootTreeWriter::GetTestVar()
    {
    	return 0.;
//    	return fTestVar;
    }

    void RootTreeWriter::SetTestVar(double aValue)
    {
//    	fTestVar = aValue;
    }

    void RootTreeWriter::OpenFile(string aFileName)
    {
        fFile = new TFile(aFileName.c_str(),"RECREATE");
    }

    void RootTreeWriter::CloseFile()
    {
    	fFile->Close();
    }


    void RootTreeWriter::WriteRootFile(Event* anEvent)
    {
    	char buffer[100];
        int n=sprintf(buffer, "Event_%d", anEvent->fEventID);
    	char* treename = buffer;

        TTree *aTree = new TTree(treename,"Locust Tree");
        aTree->Branch("EventID", &anEvent->fEventID, "EventID/I");
        aTree->Branch("ntracks", &anEvent->fNTracks, "ntracks/I");
        aTree->Branch("StartFrequencies", "std::vector<double>", &anEvent->fStartFrequencies);
        aTree->Branch("StartTimes", "std::vector<double>", &anEvent->fStartTimes);
        aTree->Branch("EndTimes", "std::vector<double>", &anEvent->fEndTimes);
        aTree->Branch("TrackLengths", "std::vector<double>", &anEvent->fTrackLengths);
        aTree->Branch("Slopes", "std::vector<double>", &anEvent->fSlopes);
        aTree->Branch("LOFrequency", &anEvent->fLOFrequency, "LOFrequency/D");
        aTree->Branch("RandomSeed", &anEvent->fRandomSeed, "RandomSeed/I");
        aTree->Branch("TrackPower", "std::vector<double>", &anEvent->fTrackPowers);
        aTree->Branch("PitchAngles", "std::vector<double>", &anEvent->fPitchAngles);
        aTree->Fill();
        aTree->Write();
        delete aTree;
    }


    //Initialize pointer to zero so that it can be initialized in first call to getInstance
//    RootTreeWriter *RootTreeWriter::instance = 0;


} /* namespace locust */
