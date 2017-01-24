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

#include "LMCGlobalsDeclaration.hh"


std::string gxml_filename = "blank.xml";

//FILE *fp2 = fopen("modeexctiation.txt","wb");  // time stamp checking.
//FILE *fp3 = fopen("fabsfakemodeexctiation.txt","wb");  // time stamp checking.


namespace locust
{
    LOGGER( lmclog, "KassSignalGenerator" );

    MT_REGISTER_GENERATOR(KassSignalGenerator, "kass-signal");

    KassSignalGenerator::KassSignalGenerator( const std::string& aName ) :
      Generator( aName ),
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
        if( aParam->Has( "xml-filename" ) )
        {
            gxml_filename = aParam->GetValue< std::string >( "xml-filename" );
        }

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
        int counter = 0;
        while (counter < 500000)  // wait 
        {
            counter += 1;
        }
    }
    if( !fKassEventReady)  // check again.
    {
        std::unique_lock< std::mutex >tLock( fMutex );
        fKassReadyCondition.wait( tLock );
        printf("LMC Got the fKassReadyCondition signal\n");
    }

    return true;

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


void* KassSignalGenerator::DriveAntenna(int PreEventCounter, unsigned index, Signal* aSignal, double* ImaginarySignal) const
{
    double c=2.998792458e8;

    locust::ParticleSlim CurrentParticle = fParticleHistory.back();
    int CurrentIndex[2]={};

    if(fParticleHistory.size()>=2 && dtConst)
    {
        double dtTmp=fParticleHistory[0].GetTime()-fParticleHistory[1].GetTime();
        if(dtConst==-99)
            dtConst=dtTmp;
        else if(fabs(dtTmp-dtConst)>1e-14)
        {
            printf("Non- constant Time Steps!");
            dtConst=0.;
        }
    }

    //Number of Grid Points Per Side. Keep odd so center point lines up with center
    const int nGridSide=9;
    double rReceiver[nGridSide][nGridSide][3];
    double rReceiverCenter[3]={0.,0.,1.};
    double ThetaSurf=0.; 
    double PhiSurf=0.;
    double tReceiverNorm[3]={sin(ThetaSurf)*cos(PhiSurf),sin(ThetaSurf)*sin(PhiSurf),cos(ThetaSurf)};
    double dx=0.005; 
    double dy=0.005;
    //Set Receiver Array Position. 
    //May eventually get from kasssiopeia 
    for(unsigned ix=0;ix<nGridSide;ix++)
    {
        for(unsigned iy=0;iy<nGridSide;iy++)
        {
            rReceiver[ix][iy][0]=double(ix)*dx-(nGridSide-1.)*dx/2.;
            rReceiver[ix][iy][1]=double(iy)*dy-(nGridSide-1.)*dy/2.;
            rReceiver[ix][iy][2]=0.;
        }
    }
    ///Perform Rotation on receiver surface
    double SurfaceRotation[3][3]={{cos(ThetaSurf)*cos(PhiSurf),-sin(PhiSurf),cos(PhiSurf)*sin(ThetaSurf)},{cos(ThetaSurf)*sin(PhiSurf),cos(PhiSurf),sin(PhiSurf)*sin(ThetaSurf)},{-sin(ThetaSurf),0.,cos(ThetaSurf)}};

    for(unsigned ix=0;ix<nGridSide;ix++)
    {
        for(unsigned iy=0;iy<nGridSide;iy++)
        {
            double tmpPos[3]={};

            for(int i=0;i<3;i++)
            {
                for(int j=0;j<3;j++)
                {
                    tmpPos[i]+=SurfaceRotation[i][j]*rReceiver[ix][iy][j];
                }
            }
            //Shift Receiver to desired location
            for(int i=0;i<3;i++)
            {
                rReceiver[ix][iy][i]=tmpPos[i]+rReceiverCenter[i];

            }
        }
    }
    

    double tRetarded=0.;
    double tReceiver=t_poststep;
    double A_Quad,B_Quad,C_Quad;
    double ReceiverPower[nGridSide][nGridSide];
    double TotalPower=0.;
    double tSpaceTimeInterval[2]={99.,99.};
    double dtRetarded[2]={};
    double tTolerance=1e-20;

    int HistorySize=int(fParticleHistory.size());
    int HistoryMaxSize=1000;
    
   

    while(tReceiver-fParticleHistory.back().GetTime()>1e-7 || fParticleHistory.size()>HistoryMaxSize)
    {
        fParticleHistory.pop_back();
        HistorySize=fParticleHistory.size();
    }
    
    printf("Size: %d\n",HistorySize);


    int AvgIters=0;

    for(unsigned ix=0;ix<nGridSide;ix++)
    {
        for(unsigned iy=0;iy<nGridSide;iy++)
        {

            if(HistorySize<HistoryMaxSize)
            {
                fParticleHistory.back().SetReceiverPosition(rReceiver[ix][iy][0],rReceiver[ix][iy][1],rReceiver[ix][iy][2]);
                fParticleHistory.back().SetReceiverTime(tReceiver);

                if(fParticleHistory.back().GetSpaceTimeInterval(0.)<0)
                {
                    ReceiverPower[ix][iy]=0.;
                    //printf("Skipping! out of Bounds!: tReceiver=%e\n",tReceiver);
                    continue;
                }
            }

            CurrentIndex[0]=FindNode(tReceiver,0.);
            CurrentParticle=fParticleHistory[CurrentIndex[0]];
            CurrentParticle.SetReceiverPosition(rReceiver[ix][iy][0],rReceiver[ix][iy][1],rReceiver[ix][iy][2]);
            CurrentParticle.SetReceiverTime(tReceiver);

            tRetarded=tReceiver-CurrentParticle.GetReceiverDistance()/c;
            //printf("tRetarded %e\n",tRetarded);

            CurrentIndex[0]=FindNode(tRetarded,CurrentIndex[0]);
            CurrentParticle=fParticleHistory[CurrentIndex[0]];
            CurrentParticle.Interpolate(tRetarded);
            CurrentParticle.SetReceiverPosition(rReceiver[ix][iy][0],rReceiver[ix][iy][1],rReceiver[ix][iy][2]);
            CurrentParticle.SetReceiverTime(tReceiver);
            tSpaceTimeInterval[0]=CurrentParticle.GetSpaceTimeInterval();

            
           // printf("Index %d of %d \n",CurrentIndex[0],HistorySize);
           // printf("Interval0: %e \n",tSpaceTimeInterval[0]);

            //CurrentParticle.CalculateQuadraticCoefficients(A_Quad,B_Quad,C_Quad);
            ////printf("A B C: %e %e %e \n",A_Quad,B_Quad,C_Quad);

            ////Roots to Quadratic Equation
            //double signB=-1.;
            //if(B_Quad>0)signB=1.;

            //dtRetarded[0]=(-B_Quad - signB *sqrt(B_Quad*B_Quad-4.*A_Quad*C_Quad)) / (2.*A_Quad);
            //dtRetarded[1]=(2.*C_Quad)/(-B_Quad- signB *sqrt(B_Quad*B_Quad-4.*A_Quad*C_Quad));

            //double OriginalIndex=CurrentIndex[0];
            //for(int i=0;i<2;i++)
            //{
            //    CurrentIndex[i]=FindNode(tRetarded+dtRetarded[i],OriginalIndex);
            //    tSpaceTimeInterval[i]=99.;
            //    //Make sure it is in bounds!
            //    if(CurrentIndex[i]>=0 && CurrentIndex[i]<HistorySize)
            //    {
            //        CurrentParticle=fParticleHistory[CurrentIndex[i]];
            //        CurrentParticle.Interpolate(tRetarded+dtRetarded[i]);
            //        CurrentParticle.SetReceiverPosition(rReceiver[ix][iy][0],rReceiver[ix][iy][1],rReceiver[ix][iy][2]);
            //        CurrentParticle.SetReceiverTime(tReceiver);
            //        tSpaceTimeInterval[i]=CurrentParticle.GetSpaceTimeInterval();
            //    }
            //}

            ////printf("dt0: %e \n",dtRetarded[0]);
            ////printf("dt1: %e \n",dtRetarded[1]);
            ////printf("s0: %e \n",tSpaceTimeInterval[0]);
            ////printf("s1: %e \n",tSpaceTimeInterval[1]);

            ////Cute way of returning index of minimum spacetime interval
            //int minIndex=int(fabs(tSpaceTimeInterval[0])>fabs(tSpaceTimeInterval[1]));
            ////printf("mindex: %d \n",minIndex);

            //dtRetarded[0]=dtRetarded[minIndex];
            //CurrentIndex[0]=CurrentIndex[minIndex];

            //tRetarded+=dtRetarded[0];
            //CurrentParticle=fParticleHistory[CurrentIndex[0]];
            //CurrentParticle.Interpolate(tRetarded);
            //CurrentParticle.SetReceiverPosition(rReceiver[ix][iy][0],rReceiver[ix][iy][1],rReceiver[ix][iy][2]);
            //CurrentParticle.SetReceiverTime(tReceiver);

            //tSpaceTimeInterval[0]=tSpaceTimeInterval[minIndex];

            //printf("s0: %e \n",tSpaceTimeInterval[0]);

            //Newton's Method

            for(int i=0;i<10;i++)
            {
                if(fabs(tSpaceTimeInterval[0])<tTolerance)
                {
                    break;
                }

                AvgIters++;

                dtRetarded[0]=tSpaceTimeInterval[0];
                //printf("s1: %e \n",dtRetarded[0]);
                //dtRetarded[0]=CurrentParticle.NewtonStep(dtRetarded[0]);
                //printf("s2: %e \n",dtRetarded[0]);
                //dtRetarded[0]=CurrentParticle.HouseHolderStep(dtRetarded[0]);
                
                tRetarded+=dtRetarded[0];
               // CurrentIndex[1]=FindNode(tRetarded,CurrentIndex[0]);

               // CurrentIndex[1]=std::max(CurrentIndex[1],0);
               // CurrentIndex[1]=std::min(CurrentIndex[1],HistorySize-1);

               // if(CurrentIndex[0]!=CurrentIndex[1])
               // {
               //     CurrentIndex[0]=CurrentIndex[1];
               //     CurrentParticle=fParticleHistory[CurrentIndex[0]];
               // }

                CurrentParticle.Interpolate(tRetarded);
                CurrentParticle.SetReceiverPosition(rReceiver[ix][iy][0],rReceiver[ix][iy][1],rReceiver[ix][iy][2]);
                CurrentParticle.SetReceiverTime(tReceiver);


                tSpaceTimeInterval[0]=CurrentParticle.GetSpaceTimeInterval();
                //printf("s1: %e \n",tSpaceTimeInterval[0]);
            }

            //printf("we got me a marlin\n");
            ReceiverPower[ix][iy]=CurrentParticle.CalculatePower(tReceiverNorm[0],tReceiverNorm[1],tReceiverNorm[2])*dx*dy;
            TotalPower+=ReceiverPower[ix][iy];

        }
    }
    printf("Power: %e\n",TotalPower);
    printf("Avg Its: %e\n",double(AvgIters)/double(nGridSide*nGridSide));


    //Scaling doesnt matter...
    double Resistance=1.;

    aLongSignal[ index ] += sqrt(TotalPower*Resistance)*cos(2.*PI*fLO_Frequency*tReceiver);
    ImaginarySignal[ index ] += sqrt(TotalPower*Resistance)*sin(2.*PI*fLO_Frequency*tReceiver);



}
////Return index of deque closest to desired time
int KassSignalGenerator::FindNode(double tNew, int IndexOld) const
{
    double tOld=fParticleHistory[IndexOld].GetTime();
    int IndexNew=0;
    bool TConst = bool(dtConst);

    if(TConst)
    {
        IndexNew=round((tOld-tNew)/dtConst)+IndexOld;
    }
    else
    {
        if(IndexOld==0)
        {
            IndexNew=BinarySearch(tNew,0,fParticleHistory.size()-1);

        }
        else
        {
            for(int i=1;i<20;i++){
                if(IsInside(tNew,IndexOld,pow(2,i)))
                {
                    IndexNew=BinarySearch(tNew,IndexOld-pow(2,i),IndexOld+pow(2,i));
                    break;
                }
            }


        }
    }

    return IndexNew;

}

bool KassSignalGenerator::IsInside(double tNew, int IndexMin, int IndexMax) const
{
    IndexMin=std::max(IndexMin,0);
    int HistorySize=int(fParticleHistory.size());
    IndexMax=std::min(IndexMax,HistorySize-1);

    return tNew>=fParticleHistory[IndexMin].GetTime() && tNew<=fParticleHistory[IndexMax].GetTime();
}

int KassSignalGenerator::BinarySearch(double tNew, int IndexMin, int IndexMax) const
{
    IndexMin=std::max(IndexMin,0);
    int HistorySize=int(fParticleHistory.size());
    IndexMax=std::min(IndexMax,HistorySize-1);

    int IndexLength=IndexMax-IndexMin;
    int IndexMid;

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
    //Now that we have it between an interval of one return index of 
    //node tNew is closer to. Ie) determine if bigger or smaller than avg. of 2 points

    return IndexMin+round(tNew-0.5*(fParticleHistory[IndexMin].GetTime()+fParticleHistory[IndexMax].GetTime()));
}




bool KassSignalGenerator::DoGenerate( Signal* aSignal ) const
{
    double *ImaginarySignal = new double[10*aSignal->TimeSize()];

    for( unsigned index = 0; index < 10*aSignal->TimeSize(); ++index )
    {
        ImaginarySignal[ index ] = 0.;  
        aLongSignal[ index ] = 0.;  // long record for oversampling.
    }

    //n samples for event spacing.
    int PreEventCounter = 0;
    int NPreEventSamples = 150000;


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
