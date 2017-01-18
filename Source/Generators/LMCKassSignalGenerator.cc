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

void* KassSignalGenerator::FilterNegativeFrequencies(Signal* aSignal) const
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
            SignalComplex[index][1] = 0.
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




void* KassSignalGenerator::DriveAntenna(int PreEventCounter, unsigned index, Signal* aSignal) const
{
    double c=2.998792458e8;
    locust::ParticleSlim CurrentParticle = fParticleHistory.front();

    //Number of Grid Points Per Side. Keep odd so center point lines up with center
    const int nGridSide=9;
    double rReceiver[nGridSide][nGridSide][3];
    double rReceiverCenter[3]={0.,0.,1.};
    double tReceiverNorm[3]={0.,0.,1.};
    double dx=0.005; 
    double dy=0.005;
    //Set Receiver Array Position. 
    //May eventually get from kasssiopeia once set experimentally
    for(unsigned ix=0;ix<nGridSide;ix++)
    {
        for(unsigned iy=0;iy<nGridSide;iy++)
        {
            rReceiver[ix][iy][0]=double(ix)*dx-(nGridSide-1.)*dx/2.;
            rReceiver[ix][iy][1]=double(iy)*dy-(nGridSide-1.)*dy/2.;
            rReceiver[ix][iy][2]=0.;
        }
    }
    //Rotate/ Shift as desired. Put in rotation in future version (just mult by matrix)
    for(unsigned ix=0;ix<nGridSide;ix++)
    {
        for(unsigned iy=0;iy<nGridSide;iy++)
        {
            for(unsigned j=0;j<3;j++)
            {
                rReceiver[ix][iy][j]+=rReceiverCenter[j];
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
    double tTolerance=1e-12;


    for(unsigned ix=0;ix<nGridSide;ix++)
    {
        for(unsigned iy=0;iy<nGridSide;iy++)
        {

            CurrentParticle.SetReceiverPosition(rReceiver[ix][iy][0],rReceiver[ix][iy][1],rReceiver[ix][iy][2]);

            //Very First Guess: Eventually Supplement with Wave Technique
            //Put into function?????
            tRetarded=tReceiver-CurrentParticle.GetReceiverDistance()/c;

            //ResetParticle(CurrentParticle,tRetarded,1);
            tSpaceTimeInterval[0]=CurrentParticle.GetSpaceTimeInterval();
            CurrentParticle.CalculateQuadraticCoefficients(A_Quad,B_Quad,C_Quad);
            //Roots to Quadratic Equation
            double signB=-1.;
            if(B_Quad>0)B_Quad=1.;

            dtRetarded[0]=(-B_Quad - signB *sqrt(B_Quad*B_Quad-4.*A_Quad*C_Quad)) / (2.*A_Quad);
            dtRetarded[1]=(2.*C_Quad)/(-B_Quad- signB *sqrt(B_Quad*B_Quad-4.*A_Quad*C_Quad));

            //Newton's Method
            while(tSpaceTimeInterval[0]>tTolerance)
            {
                CurrentParticle.CalculateQuadraticCoefficients(A_Quad,B_Quad,C_Quad);

                dtRetarded[0]=CurrentParticle.NewtonStep(tRetarded,dtRetarded[0]);
                tRetarded+=dtRetarded[0];
                CurrentIndex=FindNode(tRetarded,CurrentIndex);
            }


            ReceiverPower[ix][iy]=CurrentParticle.CalculatePower(tReceiverNorm[0],tReceiverNorm[1],tReceiverNorm[2]);
            TotalPower+=ReceiverPower[ix][iy];
            

        }
    }

    //Remove overly old elements from fParticleHistory
    while(tReceiver-fParticleHistory.back().GetTime()>3e-8 || fParticleHistory.size()>5000 )
    {
        fParticleHistory.pop_back();
    }




    //Scaling doesnt matter...
    double Resistance=1.;

    aLongSignal[ index ] += sqrt(TotalPower*Resistance);



}
////Return index of deque closest to desired time
int KassSignalGenerator::FindNode(double tNew, double tOld, int IndexOld)
{
    int IndexNew=0;
    if(TConst)
    {
        IndexNew=round((tOld-tNew)/dt)+IndexOld;
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
                if(IsInside(tNew,IndexOld,8*i))
                {
                    IndexNew=BinarySearch(tNew,IndexOld-8*i,IndexOld+8*i);
                    break;
                }
            }


        }
    }
    return IndexNew;

}

bool KassSignalGenerator::IsInside(double tNew, int IndexMin, int IndexMax)
{
    return tNew>=fParticleHistory[IndexMin].GetTime() && tNew<=fParticleHistory[IndexMax].GetTime();
}

int KassSignalGenerator::BinarySearch(double tNew, int IndexMin, int IndexMax)
{
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
    for( unsigned index = 0; index < 10*aSignal->TimeSize(); ++index )
    {
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
                    DriveAntenna(PreEventCounter, index, aSignal);
                    PreEventCounter = 0; // reset
                }
                tLock.unlock();
            }

        }  // for loop

    //FilterNegativeFrequencies(aSignal);

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
