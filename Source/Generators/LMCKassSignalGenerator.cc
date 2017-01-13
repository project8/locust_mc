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


double phi_t1 = 0.; // antenna voltage phase in radians.
double phi_t2 = 0.; // reflecting short voltage phase in radians.
double phiLO_t = 0.; // voltage phase of LO in radians;
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




//void* KassSignalGenerator::DriveAntenna(int PreEventCounter, unsigned index, Signal* aSignal) const
//{
//
//    aSignal->SignalTime()[ index ] += AverageModeExcitation()*pow(LarmorPower,0.5)*RealVoltagePhase;
//}
//

void* KassSignalGenerator::DriveAntenna(unsigned index, Signal* aSignal) const
{
    double c=2.998792458e8;

    //mass+charge
    double qParticle=qparticle;
    double mParticle=mparticle;

    //Number of Grid Points Per Side. Keep odd so center point lines up with center
    const int nGridSide=9;
    double rReceiver[nGridSide][nGridSide][3];
    double rReceiverCenter[3]={0.,0.,1.};
    double tReceiverNorm[3]={0.,0.,1.};
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

    //Solve for the retarted times for each point on receiver
    //Calculate fields directly w/ Lienard-Wiechert
    //Particle To Receiver Separation Vector/Distance
    //t are vectors
    double xParticleReceiver;
    double tParticleReceiver[3];

    double tMagneticField[3]={xMagneticField,yMagneticField,zMagneticField};
    double tMagneticFieldUnit[3]=NormalizeVector(tMagneticField);

    double rParticle[3]={X,Y,Z};

    double tParticleVelocity[3]={xVelocity,yVelocity,zVelocity};
    double tParticleVelocityUnit[3]=NormalizeVector(tParticleVelocity);
    //////////////
    double tBeta[3]={xVelocity/c,yVelocity/c,zVelocity/c};

    double xParticleVelocityParallel=Dot(tParticleVelocity,tMagneticField);
    double tParticleVelocityParallel[3]=xParticleVelocityParallel*tMagneticFieldUnit;

    double tParticleVelocityPerp[3];
    for(int i=0;i<3;i++){
        tParticleVelocityPerp[i]=tParticleVelocity[i]-tParticleVelocityParallel[i];
    }

    double xParticleVelocityPerp=Mag(tParticleVelocityPerp);

    double Gamma=1./sqrt((1.-Mag(tBeta))*(1.+Mag(tBeta)));
    double Omega=qParticle*tMagneticFieldMag/(mParticle*Gamma);
    double rCyclotron=mParticle*xParticleVelocityPerp*Gamma/(qParticle*tMagneticFieldMag);

    //Calculate Guiding Center Position:
    //May replace w/ Direct call: Right Now symplectic method doesnt have interpolation: 
    //Although could use symplectic as integrator to get as clost as possible then call with a guiding center like scheme from here
    double rParticleGuidingCenter[3];
    double tCrossRadial[3]=NormalizeVector(Cross(tParticleVelocity,tMagneticField));
    for(int i=0;i<3;i++){
        rParticleGuidingCenter[i]=rParticle[i]-rCyclotron*tCrossRadial[i];
    }

    //Other 2 vectors describing perpindicular components of motion
    double tAlpha[3], tBeta[3];
    tAlpha=rParticle-rParticleGuidingCenter;
    tBeta=Cross(tParticleVelocityUnit,tAlpha);
    NormalizeVector(tAlpha);
    NormalizeVector(tBeta);




    double tRetarded=0.;
    double tReceiver=t_poststep;
    double A_Quad,B_Quad,C_Quad;


    for(unsigned ix=0;ix<nGridSide;ix++)
    {
        for(unsigned iy=0;iy<nGridSide;iy++)
        {

            for(unsigned index=0;index<3;index++){
                tParticleReceiver[index]=rReceiver[ix][iy][index]-rParticle[index];
            }
            xParticleReceiver=Mag(tParticleReceiver);
            NormalizeVector(tParticleReceiver);

            //Very First Guess: Eventually Supplement with Wave Technique
            tRetarded=tReceiver-rParticle/c;




            

        }
    }
}

//Private Geometric Vector Functions Needed

//Cross, NormalizeVector, Mag






//double KassSignalGenerator::ModeExcitation() const
//{
//
//return EdotV;
//}
//
//
//double KassSignalGenerator::AverageModeExcitation() const
//{
//
//return AverageEdotV;
//}






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
                    DriveAntenna(PreEventCounter, index, aSignal;
                    PreEventCounter = 0; // reset
                }
                tLock.unlock();
            }

        }  // for loop


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
