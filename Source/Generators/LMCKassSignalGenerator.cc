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

#include "LMCGlobals.hh"
double phi_t = 0.; // antenna voltage phase in radians.
double phiLO_t = 0.; // voltage phase of LO;


namespace locust
{
    LOGGER( lmclog, "KassSignalGenerator" );

    MT_REGISTER_GENERATOR(KassSignalGenerator, "kass-signal");

    KassSignalGenerator::KassSignalGenerator( const std::string& aName ) :
            Generator( aName )
    {
        fRequiredSignalState = Signal::kTime;

    }

    KassSignalGenerator::~KassSignalGenerator()
    {
    }

    bool KassSignalGenerator::Configure( const ParamNode* aParam )
    {
        if( aParam == NULL) return true;

        return true;
    }

    void KassSignalGenerator::Accept( GeneratorVisitor* aVisitor ) const
    {
        aVisitor->Visit( this );
        return;
    }



    void* KassiopeiaInit()
    {
    	const string & afile = "/home/penny/project8/Project8_by_Devin.xml";
    	RunKassiopeia *RunKassiopeia1 = new RunKassiopeia;
    	RunKassiopeia1->Run(afile);
    	delete RunKassiopeia1;
    }



void WakeBeforeEvent()
{

    fPreEventCondition.notify_one();
    return;

}


bool ReceivedEventFinishSTLCondition()  // this is unused right now.
{

if( fWaitAfterEvent )
{
	printf("waiting for event to finish ...\n");
    std::unique_lock< std::mutex >tLock( fMutex );
    fPostEventCondition.wait( tLock );
    printf("event is finished in Locust\n"); getchar();
    return true;
}

return false;
}




void* DriveAntenna(unsigned index, Signal* aSignal)
{
	double dt = 5.e-9; // seconds, this should be coming from RunLengthCalculator or something like it.
    double fprime = 0.;  // Doppler shifted cyclotron frequency in Hz.

//                  printf("paused in Locust! zvelocity is %g\n", zvelocity); getchar();


           fprime = fcyc*GammaZ*(1.-zvelocity/2.99792e8);
           phi_t += 2.*PI*fprime*dt;
           phiLO_t += -2.*PI*LO_FREQUENCY*dt;

           aSignal->SignalTime()[ index ] +=
                    pow(LarmorPower,0.5)*(cos(phi_t)*cos(phiLO_t) - sin(phi_t)*sin(phiLO_t));

          printf("Locust says:  signal %d is %g and t is %g and zvelocity is %g and sqrtLarmorPower is %g and fcyc is %.10g and fprime is %g and GammaZ is %f\n",
                   index, aSignal->SignalTime()[ index ], t, zvelocity, pow(LarmorPower,0.5), fcyc, fprime, GammaZ);
//                  getchar();


}




    bool KassSignalGenerator::DoGenerate( Signal* aSignal ) const
    {

//      samples for event spacing.
    	int PreEventCounter = 0;
    	int NPreEventSamples = 15000;

    	std::thread Kassiopeia (KassiopeiaInit);     // spawn new thread
    	fRunInProgress = true;


//    	for( unsigned index = 0; index < aSignal->TimeSize(); ++index )

        for( unsigned index = 0; index < 100000; ++index )

    	{

    	if ((!fEventInProgress) && (fRunInProgress) && (!fPreEventInProgress))
    		{

    		printf("Locust says fKassEventReady is %d\n", fKassEventReady);
//    		while ((fRunInProgress) && (!fKassEventReady)) {} WakeBeforeEvent();
    		while ((fRunInProgress) && (!fKassEventReady)) {} fPreEventInProgress = true;
//    		printf("Locust just sent wakebeforeevent signal\n");
//    		if (fRunInProgress)
//    			{
  //  			fEventInProgress = true;
//    			}

    		}

    	if (fPreEventInProgress)
    	  {
    	  PreEventCounter += 1;
          printf("noise trigger index %d\n", index);
          if (PreEventCounter > NPreEventSamples)  // finished noise samples.  Start event.
        	  {
        	  fPreEventInProgress = false;  // reset.
        	  PreEventCounter = 0;  // reset.
    		  phi_t = 0.; // initialize antenna voltage phase in radians.
  		      phiLO_t = 0.; // initialize voltage phase of LO;
        	  WakeBeforeEvent();  // trigger Kass event.
              fEventInProgress = true;  // flag.
        	  }
    	  }


    	if (fEventInProgress)  // fEventInProgress
        {
//    	  printf("waiting for digitizer trigger ... index is %d\n", index);
          std::unique_lock< std::mutex >tLock( fMutexDigitizer, std::defer_lock );
          tLock.lock();
          fDigitizerCondition.wait( tLock );
          if (fEventInProgress) DriveAntenna(index, aSignal);
          tLock.unlock();
        }
    	printf("made it out of eventhappening loop. index is %d\n", index);

        }  // for loop


        // trigger any remaining events in Kassiopeia so that its thread can finish.
    	while (fRunInProgress)
    		if ((fRunInProgress) && (fKassEventReady)) WakeBeforeEvent();





//    printf("for loop is ended\n");

        Kassiopeia.join();  // wait for Kassiopeia to finish.

//   printf("kass thread has joined.\n");

        return true;
    }

} /* namespace locust */
