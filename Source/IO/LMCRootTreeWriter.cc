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


    RootTreeWriter::RootTreeWriter()
    {
    }

    RootTreeWriter::~RootTreeWriter()
    {
    }

    bool RootTreeWriter::Configure( const scarab::param_node& aParam )
    {

    	if( !FileWriter::Configure(aParam))
    	{
    		LERROR(lmclog,"Error configuring FileWriter class from RootTreeWriter child class");
    	}


    	return true;
    }

    void RootTreeWriter::WriteRunParameters( RunParameters* aRunParameter)
    {
        TTree *aTree = new TTree("Run Parameters","Locust Tree");
        aTree->Branch("LOFrequency", &aRunParameter->fLOfrequency, "LOfrequency/D");
        aTree->Branch("SamplingFrequencyMHz", &aRunParameter->fSamplingRateMHz, "SamplingFrequency/D");
        aTree->Branch("DecimationRate", &aRunParameter->fDecimationFactor, "DecimationFactor/D");
        aTree->Fill();
        aTree->Write();
        delete aTree;

    }

    void RootTreeWriter::WriteEvent(Event* anEvent)
    {
    	char buffer[100];
        int n=sprintf(buffer, "Event_%d", anEvent->fEventID);
    	char* treename = buffer;

        TTree *aTree = new TTree(treename,"Locust Tree");
        aTree->Branch("EventID", &anEvent->fEventID, "EventID/I");
        aTree->Branch("ntracks", &anEvent->fNTracks, "ntracks/I");
        aTree->Branch("StartingEnergieseV", "std::vector<double>", &anEvent->fStartingEnergies_eV);
        aTree->Branch("OutputStartFrequencies", "std::vector<double>", &anEvent->fOutputStartFrequencies);
        aTree->Branch("StartFrequencies", "std::vector<double>", &anEvent->fStartFrequencies);
        aTree->Branch("EndFrequencies", "std::vector<double>", &anEvent->fEndFrequencies);
        aTree->Branch("AvgFrequencies", "std::vector<double>", &anEvent->fAvgFrequencies);
        aTree->Branch("OutputAvgFrequencies", "std::vector<double>", &anEvent->fOutputAvgFrequencies);
        aTree->Branch("AvgAxialFrequencies", "std::vector<double>", &anEvent->fAvgAxialFrequencies);
        aTree->Branch("StartTimes", "std::vector<double>", &anEvent->fStartTimes);
        aTree->Branch("EndTimes", "std::vector<double>", &anEvent->fEndTimes);
        aTree->Branch("TrackLengths", "std::vector<double>", &anEvent->fTrackLengths);
        aTree->Branch("Slopes", "std::vector<double>", &anEvent->fSlopes);
        aTree->Branch("LOFrequency", &anEvent->fLOFrequency, "LOFrequency/D");
        aTree->Branch("RandomSeed", &anEvent->fRandomSeed, "RandomSeed/I");
        aTree->Branch("TrackPower", "std::vector<double>", &anEvent->fTrackPowers);
        aTree->Branch("PitchAngles", "std::vector<double>", &anEvent->fPitchAngles);
        aTree->Branch("Radii", "std::vector<double>", &anEvent->fRadii);
        aTree->Branch("RadialPhases", "std::vector<double>", &anEvent->fRadialPhases);
        aTree->Fill();
        aTree->Write();
        delete aTree;
    }



} /* namespace locust */
