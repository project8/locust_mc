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
        TransmitterInterfaceGenerator( aName )
    {
    }

    FieldParameterGenerator::~FieldParameterGenerator()
    {
    }

    bool FieldParameterGenerator::Configure( const scarab::param_node& aParam )
    {
	    TransmitterInterfaceGenerator::Configure(aParam);
    }

    void FieldParameterGenerator::InitializeFieldPoints()
    {
	for(int pointIndex = 0; pointIndex< GetNPoints(); ++pointIndex)
	{
	    LMCThreeVector point(0.0,0.0,0.0);
	    InitializeFieldPoint(point);
	}
    }

    static void* KassiopeiaInit(const std::string &aFile)
    {
        RunKassiopeia RunKassiopeia1;
        RunKassiopeia1.Run(aFile);
        RunKassiopeia1.~RunKassiopeia();

        return 0;
    }

    void FieldParameterGenerator::DriveAntenna(FILE *fp, int PreEventCounter, unsigned index, Signal* aSignal, int nfilterbins, double dtfilter)
    {
        const int signalSize = aSignal->TimeSize();
        unsigned pointIndex = 0;
        unsigned sampleIndex = 0;

        unsigned tTotalElementIndex = 0;

        for(int pointIndex = 0; pointIndex < GetNPoints(); ++pointIndex)
        {
                sampleIndex = pointIndex*aSignal->DecimationFactor() + index;  // which point and which sample

                double* tFieldSolution = new double[2];
                if (!fTransmitter->IsKassiopeia())
                {
                	tFieldSolution = fTransmitter->GetEFieldCoPol(pointIndex, 1./(fAcquisitionRate*1.e6*aSignal->DecimationFactor()));
                }
                else
                {
                	tFieldSolution = fTransmitter->SolveKassFields(fAllFieldCopol[pointIndex],fAllFieldCopol[pointIndex],t_old,pointIndex);
                }
                if (fTextFileWriting==1) {}
		RecordIncidentFields(fp, t_old,fAllFieldCopol.at(pointIndex), tFieldSolution[1]);
 	        FillBuffers(aSignal, tFieldSolution[1], tFieldSolution[0],pointIndex,index);
                PopBuffers(pointIndex);
        } // channels loop

        t_old += 1./(fAcquisitionRate*1.e6*aSignal->DecimationFactor());
        if ( index%fSwapFrequency == 0 ) CleanupBuffers();  // release memory
    }

    bool FieldParameterGenerator::DoGenerate( Signal* aSignal )
    {
        FILE *fp = fopen("incidentfields.txt", "w");

        //n samples for event spacing in Kass.
        int PreEventCounter = 0;
        const int NPreEventSamples = 150000;
        fKassTimeStep = 1./(fAcquisitionRate*1.e6*aSignal->DecimationFactor());

        InitializeBuffers(fNFilterBins, fFieldBufferSize);
        InitializeFieldPoints();

        if (!fTransmitter->IsKassiopeia())
        {
        	for( unsigned index = 0; index < aSignal->DecimationFactor()*aSignal->TimeSize(); ++index )
        	{
        		DriveAntenna(fp, PreEventCounter, index, aSignal, fNFilterBins,fdtFilter);
        	}  // for loop
        	return true;
        }

        else if (fTransmitter->IsKassiopeia())
        {

            std::thread Kassiopeia(KassiopeiaInit, gxml_filename);  // spawn new thread
        	fRunInProgress = true;

        for( unsigned index = 0; index < aSignal->DecimationFactor()*aSignal->TimeSize(); ++index )
        {
            if ((!fEventInProgress) && (fRunInProgress) && (!fPreEventInProgress))
            {
                if (ReceivedKassReady()) fPreEventInProgress = true;
            }

            if (fPreEventInProgress)
            {
                PreEventCounter += 1;
                //printf("preeventcounter is %d\n", PreEventCounter);
                if (PreEventCounter > NPreEventSamples)  // finished noise samples.  Start event.
                {
                    fPreEventInProgress = false;  // reset.
                    fEventInProgress = true;
                    //printf("LMC about to wakebeforeevent\n");
                    WakeBeforeEvent();  // trigger Kass event.
                }
            }

            if (fEventInProgress)  // fEventInProgress
                if (fEventInProgress)  // check again.
                {
                    //printf("waiting for digitizer trigger ... index is %d\n", index);
                    std::unique_lock< std::mutex >tLock( fMutexDigitizer, std::defer_lock );
                    tLock.lock();
                    fDigitizerCondition.wait( tLock );
                    if (fEventInProgress)
                    {
                        //printf("about to drive antenna, PEV is %d\n", PreEventCounter);
                        DriveAntenna(fp, PreEventCounter, index, aSignal, fNFilterBins,fdtFilter);
                        PreEventCounter = 0; // reset
                    }
                    tLock.unlock();
                }

        }  // for loop

        printf("finished signal loop\n");

        fclose(fp);
        fRunInProgress = false;  // tell Kassiopeia to finish.
        fDoneWithSignalGeneration = true;  // tell LMCCyclotronRadExtractor
        WakeBeforeEvent();
        Kassiopeia.join();

        return true;
        }

    }

} /* namespace locust */

