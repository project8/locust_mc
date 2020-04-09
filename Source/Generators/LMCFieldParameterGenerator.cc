/*
 * LMCFieldParameterGenerator.cc
 *
 *  Created on: Mar 04, 2020
 *      Author: P. T. Surukuchi
 */

#include "LMCFieldParameterGenerator.hh"
#include "LMCSphere.hh"

namespace locust
{
    LOGGER( lmclog, "FieldParameterGenerator" );

    MT_REGISTER_GENERATOR(FieldParameterGenerator, "field-parameter");

    FieldParameterGenerator::FieldParameterGenerator( const std::string& aName ):
        TransmitterInterfaceGenerator( aName ),
        fUseTextFile(false),
        fPredefinedGeometry(0),
        fRadius(1.0),
        fLength(0.0),
        fMinFieldPoints(0),
        fTextFileName("blank.txt")
    {
    }

    FieldParameterGenerator::~FieldParameterGenerator()
    {
    }

    bool FieldParameterGenerator::ends_with(const std::string &str, const std::string &suffix)
    {
        //copied from https://stackoverflow.com/a/20446239
        return str.size() >= suffix.size() &&
            str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
    }

    bool FieldParameterGenerator::Configure( const scarab::param_node& aParam )
    {
        TransmitterInterfaceGenerator::Configure(aParam);
        if(aParam.has("field-point-inputs"))
        {
            const scarab::param_node& aParamFieldPoints=aParam["field-point-inputs"].as_node();
            if(aParamFieldPoints.has("use-text-file"))
            {
                fUseTextFile=aParamFieldPoints["use-text-file"]().as_bool();
                fTextFileName=aParamFieldPoints["file-name"]().as_string();
            }
            if(fUseTextFile)
            {
                fUseTextFile=true;
                if(aParamFieldPoints.has("file-name"))
                {
                    LDEBUG(lmclog,"Using the text file "<<fTextFileName<< " for defining field points in the FieldParameterGenerator");
                }
                else
                {
                    LERROR(lmclog,"A filename with the config option 'file-name' has to be input through the config file to define field points");
                    exit(-1);
                }
            }
            else
            {
                std::string aPredefinedGeometry="sphere";
                if(aParamFieldPoints.has("predefined-geometry"))
                {
                    aPredefinedGeometry=aParamFieldPoints["predefined-geometry"]().as_string();
                    if(aParamFieldPoints.has("radius"))
                    {
                        fRadius=aParamFieldPoints["radius"]().as_double();
                    }
                    if(aParamFieldPoints.has("minimum-field-points"))
                    {
                        fMinFieldPoints=aParamFieldPoints["minimum-field-points"]().as_double();
                    }
                    if(aPredefinedGeometry.compare("sphere")==0) fPredefinedGeometry=0;
                    else if(aPredefinedGeometry.compare("cylinder")==0) 
                    {
                        fPredefinedGeometry=1;
                        if(aParamFieldPoints.has("length"))
                        {
                            fLength=aParamFieldPoints["length"]().as_double();
                        }
                    }
                    else
                    {
                        LERROR(lmclog,"The geometry type for "<< aPredefinedGeometry<< " is not currenlty defined in FieldParameterGenerator");
                        exit(-1);
                    }
                }
                LDEBUG(lmclog,"Using a predefined geometry of "<<aPredefinedGeometry.c_str()<< " for defining field points in the FieldParameterGenerator");
            }
        }
        else
        {
            fPredefinedGeometry=0;
            LWARN(lmclog,"Using a predefined spherical geometry of radius 1 for defining field points in the FieldParameterGenerator");
        }
        InitializeFieldPoints();
        return true;
    }

    double FieldParameterGenerator::GenerateFieldPoints()
    {
        if(fUseTextFile)
        {
            if(!ends_with(fTextFileName.c_str(),".txt"))
            {
                LERROR(lmclog,"The text file " << fTextFileName.c_str() <<" doesn't end in .txt");
                exit(-1);
            }
            double xPos;
            double yPos;
            double zPos;
            std::fstream textFile(fTextFileName.c_str(),std::ios::in);
            if (textFile.fail())
            {
                LERROR(lmclog,"The file " <<fTextFileName.c_str() <<" doesn't exist");
                exit(-1);
            }
            while(!textFile.eof())
            {
                std::string lineContent;
                while(std::getline(textFile,lineContent))
                {
                    if (lineContent.find('#')!=std::string::npos) continue;
                    std::string token;
                    std::stringstream ss(lineContent);
                    int wordCount=0;
                    while (ss >> token)
                    {
                        if(wordCount==0)xPos=std::stod(token);
                        else if(wordCount==1)yPos=std::stod(token);
                        else if(wordCount==2)zPos=std::stod(token);
                        else
                        {
                            LERROR(lmclog, "The input file for field points should only have three columns");
                            exit(-1);
                        }
                        ++wordCount;
                    }
                    fFieldPoints.push_back(LMCThreeVector(xPos,yPos,zPos));
                }
            }
        }
        else 
        {
            if(fPredefinedGeometry==0)
            {
                LMCIcoSphere sphere(fRadius,LMCThreeVector(),fMinFieldPoints);
                sphere.GetVertices(fFieldPoints); 
            }
            else if(fPredefinedGeometry==1)
            {
                LERROR(lmclog,"FieldParameterGenerator doesn't have cylinder implementation yet");
            }
        }
        SetNPoints(fFieldPoints.size());
        LDEBUG(lmclog,"FieldParameterGenerator built with "<<GetNPoints()<< " points");
        return GetNPoints();
    }

    void FieldParameterGenerator::InitializeFieldPoints()
    {
        GenerateFieldPoints();
        for(int pointIndex = 0; pointIndex< GetNPoints(); ++pointIndex)
        {
            InitializeFieldPoint(fFieldPoints.at(pointIndex));
        }
        exit(-1);
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
            std::cout<<tFieldSolution[1]<<"  "<< tFieldSolution[0]<< "  "<<pointIndex<< " "<< index<<std::endl;
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

