/*
 * LMCKassSignalGenerator.cc
 *
 *  Created on: Mar 12, 2014
 *      Author: nsoblath
 */

#include "LMCKassSignalGenerator.hh"
#include "LMCEventHold.hh"
#include "LMCRunKassiopeia.hh"

#include "logger.hh"
#include <thread>
#include <algorithm>

#include <iostream>
#include <fstream>

#include "LMCGlobalsDeclaration.hh"
#include "LMCHFSSReader.hh"


std::string gxml_filename = "blank.xml";

double phi_LO=0.;
const int NPreEventSamples = 150000;
int fNFDIndex=-1;

//double fPowerZ[100]={};

namespace locust
{
    LOGGER( lmclog, "KassSignalGenerator" );

    MT_REGISTER_GENERATOR(KassSignalGenerator, "kass-signal");

    KassSignalGenerator::KassSignalGenerator( const std::string& aName ) :
      Generator( aName ),
      fWriteNFD(0.),
      fLO_Frequency( 0.)
    {
        fRequiredSignalState = Signal::kTime;
    }

    KassSignalGenerator::~KassSignalGenerator()
    {
    }

    bool KassSignalGenerator::Configure( const ParamNode* aParam )
    {
        if( aParam == NULL) return true;
        fLO_Frequency = LO_FREQUENCY;
        if( aParam->Has( "lo-frequency" ) )
        {
            fLO_Frequency = aParam->GetValue< double >( "lo-frequency" );
        }
        if( aParam->Has( "carrier-frequency" ) )
        {
            double fCarrier_Frequency = aParam->GetValue< double >( "carrier-frequency" );
            fDecimationFactor = AntiAliasingSetup(fCarrier_Frequency,1.e9);
            DEBUG( lmclog, "Changing Decimation Factor to: "<< fDecimationFactor);
            EventModTimeStep*=10./double(fDecimationFactor);
        }
        if( aParam->Has( "xml-filename" ) )
        {
            gxml_filename = aParam->GetValue< std::string >( "xml-filename" );
        }
        if( aParam->Has( "and-filename" ) )
        {
            fWriteNFD=1;
            fAND_filename = aParam->GetValue< std::string >( "and-filename" );
            HFSSReader HFRead;
            HFRead.ParseANDFile(fAND_filename);
            rReceiver=HFRead.GetSurfacePoints();
            NFDFrequencies=HFRead.GetFrequencies();
            fNFD_filename=HFRead.GetNFDFilename();
        }
        else
        {
            //If not, use existing code to generate plane receiver
            HFSSReader HFRead;
            rReceiver=HFRead.GeneratePlane({0.05,0.05},7);//Argumemts: Size, resolution
            rReceiver=HFRead.RotateShift(rReceiver,{1.,0.,0.},{0.05,0.,0.});//Arguments Normal vector, Position (m)
        }
        PreviousTimes=std::vector<std::array<double,2> >(rReceiver.size(),{-99.,-99.});

        return true;
    }

    void KassSignalGenerator::Accept( GeneratorVisitor* aVisitor ) const
    {
        aVisitor->Visit( this );
        return;
    }


  void* KassiopeiaInit()
    {
        //cout << gxml_filename; getchar();
        const string & afile = gxml_filename;
    	RunKassiopeia *RunKassiopeia1 = new RunKassiopeia;
    	RunKassiopeia1->Run(afile);
    	delete RunKassiopeia1;
    }



void WakeBeforeEvent()
{

    fPreEventCondition.notify_one();
    return;

}




bool ReceivedKassReady()
{

    if( !fKassEventReady)
    {
    std::unique_lock< std::mutex >tLock( fKassReadyMutex );
    fKassReadyCondition.wait( tLock );
    printf("LMC Got the fKassReadyCondition signal\n");
    }

    return true;

}



//Change decimation factor/ (and therefore the sampling rate) to guarantee no aliasing of signal
double KassSignalGenerator::AntiAliasingSetup(double fCarrier_Frequency, double fBandwidth_Frequency) const
{
    double fSampleMin = 2.*fBandwidth_Frequency;
    int fDecimationRange[2] = {10,100};
    int fDecimation = fDecimationRange[0];

    double fSample_Frequency;
    int mFactor;

    while(fDecimation <= fDecimationRange[1])
    {
        fSample_Frequency=double(fDecimation)/double(fDecimationRange[0])*fSampleMin;
        mFactor=floor((2.*fCarrier_Frequency-fBandwidth_Frequency)/fSample_Frequency);
        if(fSample_Frequency >= (2.*fCarrier_Frequency+fBandwidth_Frequency)/(double(mFactor+1)))
        {
            break;
        }

        fDecimation++;
    }

    if(fDecimation==fDecimationRange[1])
    {
            ERROR( lmclog, "Cannot find Decimation Factor. Are you sure about your carrier frequency?");
    }

    return fDecimation;
}

void* KassSignalGenerator::FilterNegativeFrequencies(Signal* aSignal, double *ImaginarySignal) const
{

    int nwindows = 80;
    int windowsize = 10*aSignal->TimeSize()/nwindows;

	fftw_complex *SignalComplex;
    SignalComplex = (fftw_complex*)fftw_malloc( sizeof(fftw_complex) * windowsize );

	fftw_complex *FFTComplex;
    FFTComplex = (fftw_complex*)fftw_malloc( sizeof(fftw_complex) * windowsize );

    fftw_plan ForwardPlan;
    ForwardPlan = fftw_plan_dft_1d(windowsize, SignalComplex, FFTComplex, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan ReversePlan;
    ReversePlan = fftw_plan_dft_1d(windowsize, FFTComplex, SignalComplex, FFTW_BACKWARD, FFTW_ESTIMATE);

    for (int nwin = 0; nwin < nwindows; nwin++)
    {
        //Put (Real Voltage into Time Domain)
        for( unsigned index = 0; index < windowsize; ++index )
        {
            SignalComplex[index][0] = aLongSignal[ nwin*windowsize + index ];
            SignalComplex[index][1] = ImaginarySignal[ nwin*windowsize + index ];
        }

        fftw_execute(ForwardPlan);
        //Complex filter to set power at negative frequencies to 0.
        for( unsigned index = windowsize/2; index < windowsize; ++index )
        {
            FFTComplex[index][0] = 0.;
            FFTComplex[index][1] = 0.;
        }

        fftw_execute(ReversePlan);

        double norm = (double)(windowsize);

        for( unsigned index = 0; index < windowsize; ++index )
        {
            //normalize and take the real part of the reverse transform, for digitization.
            aLongSignal[ nwin*windowsize + index ] = SignalComplex[index][0]/norm;
        }

    }


    delete SignalComplex;
    delete FFTComplex;

}

void KassSignalGenerator::NFDWrite() const
{
        std::ofstream fNFDOutput;
        fNFDOutput.open(fNFD_filename,std::ios::out | std::ios::trunc);

        fNFDOutput << "#Index, X, Y, Z, Ex(real, imag), Ey(real, imag), Ez(real, imag), Hx(real, imag), Hy(real, imag), Hz(real, imag)"<<std::endl;
        
        fNFDOutput << std::fixed;
        fNFDOutput << "Frequencies "<<NFDFrequencies.size()<<std::endl;

        fNFDOutput << std::setprecision(9);

        // Loop over frequencies
        for(int i = 0; i < NFDFrequencies.size() ; i++)
        {

            fNFDOutput<<std::scientific;
            fNFDOutput <<"Frequency "<< NFDFrequencies[i] <<std::endl;

            for(int j = 0; j < rReceiver.size() ; j++)
            {

                fNFDOutput<<std::scientific;
                fNFDOutput <<j+1<<", "<<rReceiver[j][0]<<", "<<rReceiver[j][1]<<", "<<rReceiver[j][2]<<", "<< NFDElectricFieldFreq[i][j][0][0]<<", "<< NFDElectricFieldFreq[i][j][0][1]<<", "<< NFDElectricFieldFreq[i][j][1][0]<<", "<< NFDElectricFieldFreq[i][j][1][1]<<", "<< NFDElectricFieldFreq[i][j][2][0]<<", "<< NFDElectricFieldFreq[i][j][2][1]<<", "<< NFDMagneticFieldFreq[i][j][0][0]<<", "<< NFDMagneticFieldFreq[i][j][0][1]<<", "<< NFDMagneticFieldFreq[i][j][1][0]<<", "<< NFDMagneticFieldFreq[i][j][1][1]<<", "<< NFDMagneticFieldFreq[i][j][2][0]<<", "<< NFDMagneticFieldFreq[i][j][2][1]<<std::endl;

            }
        }
        fNFDOutput.close();


        //Be 100% sure vector deallocates memory correctly (probably overkill)
        NFDElectricFieldFreq.resize(0);
        NFDElectricFieldFreq.shrink_to_fit();

        NFDMagneticFieldFreq.resize(0);
        NFDMagneticFieldFreq.shrink_to_fit();

}


void* KassSignalGenerator::DriveAntenna(int PreEventCounter, unsigned index, Signal* aSignal, double* ImaginarySignal) const
{
    //std::ofstream myfile;
    //myfile.open ("Poynting.txt",std::ios::app);

    locust::ParticleSlim CurrentParticle = fParticleHistory.back();
    int CurrentIndex;

    double tReceiver=t_old;
    double tRetarded=0.;
    //double tReceiver=t_poststep;
    double ReceiverVoltage[rReceiver.size()];
    double TotalVoltage=0.;

    double tSpaceTimeInterval=99.;
    double dtRetarded=0;
    double tTolerance=1e-23;
    double dtStepSize=fabs(fParticleHistory[1].GetTime()-fParticleHistory[0].GetTime());

    int HistorySize=fParticleHistory.size();
    int HistoryMaxSize=5000;

    phi_LO+=2.*PI*fLO_Frequency*EventModTimeStep;

   //printf("Size: %d %d\n",HistorySize, fParticleHistory.size());

    int AvgIters=0;

    for(unsigned i=0;i<rReceiver.size();i++)
    {
        //Check if there is time for photon to reach receiver if particle is recently created
        if(fParticleHistory.front().GetTime()<=3.*dtStepSize)
        {
            fParticleHistory.front().SetReceiverPosition(rReceiver[i][0],rReceiver[i][1],rReceiver[i][2]);
            fParticleHistory.front().SetReceiverTime(tReceiver);

            if(fParticleHistory.front().GetSpaceTimeInterval(0.)<0)
            {
                ReceiverVoltage[i]=0.;
                printf("Skipping! out of Bounds!: tReceiver=%e\n",tReceiver);
                continue;
            }
        }

        if(PreviousTimes[i][0]==-99.)
        {
            CurrentIndex=FindNode(tReceiver,dtStepSize,HistorySize-1);
            CurrentParticle=fParticleHistory[CurrentIndex];
            CurrentParticle.SetReceiverPosition(rReceiver[i][0],rReceiver[i][1],rReceiver[i][2]);
            CurrentParticle.SetReceiverTime(tReceiver);

            double c=2.998792458e8;
            tRetarded=tReceiver-CurrentParticle.GetReceiverDistance()/c;

            //printf("Index %d of %d\n",CurrentIndex[0],HistorySize-1);
            //printf("--------\n");
        }
        else
        {
            //printf("EventModStep: %e\n",EventModTimeStep);
            tRetarded=PreviousTimes[i][0]+EventModTimeStep;
            CurrentIndex=int(PreviousTimes[i][1]);
        }

        CurrentIndex=std::max(CurrentIndex,0);
        CurrentIndex=std::min(CurrentIndex,HistorySize-1);

        CurrentIndex=FindNode(tRetarded,dtStepSize,CurrentIndex);

        CurrentParticle=fParticleHistory[CurrentIndex];
        CurrentParticle.Interpolate(tRetarded);
        CurrentParticle.SetReceiverPosition(rReceiver[i][0],rReceiver[i][1],rReceiver[i][2]);
        CurrentParticle.SetReceiverTime(tReceiver);
        tSpaceTimeInterval=CurrentParticle.GetSpaceTimeInterval();

        //printf("s0: %e \n",tSpaceTimeInterval);

        double tOldSpaceTimeInterval=99.;

        //Converge to root
        for(int j=0;j<25;j++)
        {

            AvgIters++;

            tRetarded=CurrentParticle.GetStepRoot(0,tSpaceTimeInterval);

            CurrentParticle.Interpolate(tRetarded);

            if(fabs(CurrentParticle.GetTimeDisplacement()) > dtStepSize)
            {
                CurrentIndex=FindNode(tRetarded,dtStepSize,CurrentIndex);
                CurrentParticle=fParticleHistory[CurrentIndex];
                CurrentParticle.Interpolate(tRetarded);
            }
            //printf("%e %e\n",CurrentParticle.GetTimeDisplacement(),dtStepSize*0.75);

            CurrentParticle.SetReceiverPosition(rReceiver[i][0],rReceiver[i][1],rReceiver[i][2]);
            CurrentParticle.SetReceiverTime(tReceiver);

            tSpaceTimeInterval=CurrentParticle.GetSpaceTimeInterval();

            tOldSpaceTimeInterval=tSpaceTimeInterval;
        }

        PreviousTimes[i][0]=tRetarded;
        PreviousTimes[i][1]=double(CurrentIndex);

        ReceiverVoltage[i]=CurrentParticle.CalculateVoltage();
        TotalVoltage+=ReceiverVoltage[i];

        double tMinHFSS=1e-8;
        const int nHFSSBins=2048;
        if(fWriteNFD && fNFDIndex<nHFSSBins && tReceiver >= tMinHFSS )
        {
            double tmpElectricField[3]; double tmpMagneticField[3];
            CurrentParticle.CalculateElectricField(tmpElectricField[0],tmpElectricField[1],tmpElectricField[2]);
            CurrentParticle.CalculateMagneticField(tmpMagneticField[0],tmpMagneticField[1],tmpMagneticField[2]);
            
            //If doing the first receiver
            if(i == 0) fNFDIndex++;


            for(int j=0;j<NFDFrequencies.size();j++)
            {
                for(int k=0;k<3;k++)
                {
                    //Complex factors for Downmixing/ DFT sums
                    double DownConvert[2]={cos(phi_LO),-sin(phi_LO)};
                    int kFreq=int((NFDFrequencies[j]-fLO_Frequency)*double(nHFSSBins)*EventModTimeStep);
                    double DFTFactor[2]={cos(2.*PI*double(kFreq)/double(nHFSSBins)*double(fNFDIndex)),-sin(2.*PI*double(kFreq)/double(nHFSSBins)*double(fNFDIndex))};

                    NFDElectricFieldFreq[j][i][k][0]+=tmpElectricField[k]*(DownConvert[0]*DFTFactor[0]-DownConvert[1]*DFTFactor[1])/(double(nHFSSBins));
                    NFDElectricFieldFreq[j][i][k][1]+=tmpElectricField[k]*(DownConvert[0]*DFTFactor[1]+DownConvert[1]*DFTFactor[0])/(double(nHFSSBins));

                    NFDMagneticFieldFreq[j][i][k][0]+=tmpMagneticField[k]*(DownConvert[0]*DFTFactor[0]-DownConvert[1]*DFTFactor[1])/(double(nHFSSBins));
                    NFDMagneticFieldFreq[j][i][k][1]+=tmpMagneticField[k]*(DownConvert[0]*DFTFactor[1]+DownConvert[1]*DFTFactor[0])/(double(nHFSSBins));
                }
            }
            if(fNFDIndex%100==0 && i==0)printf("%d\n",fNFDIndex);

        }
        else if(fWriteNFD && fNFDIndex >=nHFSSBins)
        {
            printf("Done! \n");
        }
    }

    double Voltage=TotalVoltage/double(rReceiver.size());

    aLongSignal[ index ] += Voltage*cos(phi_LO);
    ImaginarySignal[ index ] += -Voltage*sin(phi_LO);

}
////Return index of deque closest to desired time
int KassSignalGenerator::FindNode(double tNew, double dtStepSize, int IndexOld) const
{
    int HistorySize=fParticleHistory.size();

    IndexOld=std::max(IndexOld,0);
    IndexOld=std::min(IndexOld,HistorySize-1);

    double tOld=fParticleHistory[IndexOld].GetTime();
    int IndexMid=round((tNew-tOld)/dtStepSize)+IndexOld;

    IndexMid=std::max(IndexMid,0);
    IndexMid=std::min(IndexMid,HistorySize-1);

    double tMid=fParticleHistory[IndexMid].GetTime();

    int IndexNew=0;

    for(int i=0;i<20;i++){
        if(IsInside(tNew,IndexMid-pow(2,i),IndexMid+pow(2,i)))
        {
            IndexNew=BinarySearch(tNew,IndexMid-pow(2,i),IndexMid+pow(2,i));
            //printf("New Index: %d\n",IndexNew);
            break;
        }
    }

    return IndexNew;

}

bool KassSignalGenerator::IsInside(double tNew, int IndexMin, int IndexMax) const
{
    int HistorySize=fParticleHistory.size();
    IndexMin=std::max(IndexMin,0);
    IndexMax=std::min(IndexMax,HistorySize-1);

    return tNew>=fParticleHistory[IndexMin].GetTime() && tNew<=fParticleHistory[IndexMax].GetTime();
}

int KassSignalGenerator::BinarySearch(double tNew, int IndexMin, int IndexMax) const
{
    int HistorySize=fParticleHistory.size();
    IndexMin=std::max(IndexMin,0);
    IndexMax=std::min(IndexMax,HistorySize-1);

    int IndexLength=IndexMax-IndexMin;
    int IndexMid=0;

    while(IndexLength>1)
    {
        IndexMid=round((IndexMin+IndexMax)/2.);
        if(fParticleHistory[IndexMid].GetTime()>tNew)
        {
            IndexMax=IndexMid;

        }
        else if(fParticleHistory[IndexMid].GetTime()<tNew)
        {
            IndexMin=IndexMid;
        }
        else if(fParticleHistory[IndexMid].GetTime()==tNew)
        {
            return IndexMid;
        }

        IndexLength=IndexMax-IndexMin;
    }
    
    IndexMin=std::max(IndexMin,0);
    return IndexMin;
}

bool KassSignalGenerator::DoGenerate( Signal* aSignal ) const
{
    double *ImaginarySignal = new double[10*aSignal->TimeSize()];

    for( unsigned index = 0; index < 10*aSignal->TimeSize(); ++index )
    {
        ImaginarySignal[ index ] = 0.;
        aLongSignal[ index ] = 0.;  // long record for oversampling.
    }

    if(fWriteNFD)
    {
        //Initialize to correct sizes/ all zeros
        NFDElectricFieldFreq.resize(NFDFrequencies.size());
        NFDMagneticFieldFreq.resize(NFDFrequencies.size());
        for(int i=0;i<NFDFrequencies.size();i++)
        {
            NFDElectricFieldFreq[i]=std::vector<std::array<std::array<double,2>, 3> >(rReceiver.size(),{0.});
            NFDMagneticFieldFreq[i]=std::vector<std::array<std::array<double,2>, 3> >(rReceiver.size(),{0.});
        }
        
    }

    //n samples for event spacing.
    int PreEventCounter = 0;

    std::thread Kassiopeia (KassiopeiaInit);     // spawn new thread
    fRunInProgress = true;

	for( unsigned index = 0; index < 10*aSignal->TimeSize(); ++index )
    {
        if ((!fEventInProgress) && (fRunInProgress) && (!fPreEventInProgress))
        	{
                if (ReceivedKassReady()) fPreEventInProgress = true;
        	}

    	if (fPreEventInProgress)
        {
            PreEventCounter += 1;
            // printf("preeventcounter is %d\n", PreEventCounter);
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
                    DriveAntenna(PreEventCounter, index, aSignal, ImaginarySignal);

                    PreEventCounter = 0; // reset
                }
                tLock.unlock();
            }

        }  // for loop

    if(fWriteNFD) NFDWrite();
    
    FilterNegativeFrequencies(aSignal, ImaginarySignal);
    delete ImaginarySignal;
    
    // trigger any remaining events in Kassiopeia so that its thread can finish.
    while (fRunInProgress)
    {
        if (fRunInProgress)
            if (ReceivedKassReady()) WakeBeforeEvent();
    }

    Kassiopeia.join();  // wait for Kassiopeia to finish.

    return true;
    }

} /* namespace locust */
