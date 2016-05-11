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




void* KassSignalGenerator::DriveAntenna(unsigned index, Signal* aSignal) const
{
	double dummy = 0.;
	double dt = 5.e-9; // seconds, this should be coming from RunLengthCalculator or something like it.
    double fprime = 0.;  // Doppler shifted cyclotron frequency in Hz.

//                  printf("paused in Locust! zvelocity is %g\n", zvelocity); getchar();


           fprime = fcyc*GammaZ*(1.-zvelocity/2.99792e8);
           phi_t += 2.*PI*fprime*dt;
           phiLO_t += -2.*PI*LO_FREQUENCY*dt;

           aSignal->SignalTime()[ index ] +=
                    ModeExcitation()*pow(LarmorPower,0.5)*(cos(phi_t)*cos(phiLO_t) - sin(phi_t)*sin(phiLO_t));

           printf("driving antenna, ModeExcitation is %g\n", ModeExcitation());
           printf("Locust says:  signal %d is %g and t is %g and zvelocity is %g and sqrtLarmorPower is %g and fcyc is %.10g and fprime is %g and GammaZ is %f\n",
                   index, aSignal->SignalTime()[ index ], t_poststep, zvelocity, pow(LarmorPower,0.5), fcyc, fprime, GammaZ);
//                  getchar();

}

double KassSignalGenerator::ModeExcitation() const
{
    double dim1_wr42 = 10.668e-3; // m
    double dim2_wr42 = 4.318e-3; // m
	double vy = yvelocity;
	double vx = xvelocity;
	double x = X + dim1_wr42/2.;  // center of waveguide is at zero.
	double Ey = 0.;
	double EyMax = 0.;

	double *EyArray1 = EyWR42Array();

//  normalize Ey if necessary.
//	EyArray1 = ScaleArray(EyArray1, 1./pow(IntEyWR42ArraySqdA(EyArray1, dim1_wr42, dim2_wr42), 0.5));

//	x=0.+dim1_wr42/2.; vy = 5.e7; vx=0.;  // fake test calcs in middle of waveguide.


	Ey = EyArray1[(int)(100.*x/dim1_wr42)];
	EyMax = EyArray1[100/2];
//	printf("EyMax is %g\n", EyMax);



	// E dot v / (Emax v) / sqrt(2) for half power lost in opposite direction.
	double EdotV = fabs(Ey*vy)/fabs(EyMax*pow(vx*vx+vy*vy,0.5)) / 1.41421;
	printf("x is %f and Ey is %g and yvelocity is %g and xvelocity is %g and EdotV is %f\n", x, Ey, vy, vx, EdotV);

return EdotV;
}


double* KassSignalGenerator::EyWR42Array() const
{
double a = 10.668e-3;
int nbins = 100;
double *EyArray1 = new double[nbins];
double x=0.;
for (int i=0; i<nbins; i++)
  {
  x = a*(double)i/(double)nbins;
  EyArray1[i] = sin(PI*x/a);
  }
return EyArray1;
}


double KassSignalGenerator::IntEyWR42ArraySqdA(double *EyArray1, double dimx, double dimy) const
{
double Int = 0.;
int nbins = 100;
double dx = dimx/(double)nbins;

for (int i=0; i<nbins; i++)
  {
  Int += EyArray1[i]*EyArray1[i]*dx;
//  printf("int is %f\n", Int);
  }

//  printf("int of histo is %f\n", Int);

Int *= dimy;
return Int;
}


double* KassSignalGenerator::ScaleArray(double *array, double factor) const
{
int nbins = 100;
for (int i=0; i<nbins; i++)
  array[i] *= factor;
return array;
}




    bool KassSignalGenerator::DoGenerate( Signal* aSignal ) const
    {

    	testvar = 5.;


//      samples for event spacing.
    	int PreEventCounter = 0;
    	int NPreEventSamples = 15000;

    	FILE *fp = fopen("timing.txt","wb");
//    	fprintf(fp, "testing\n");


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
          if (fEventInProgress)
        	  {
        	  fprintf(fp, "%.10g  ", t_poststep);
        	  DriveAntenna(index, aSignal);
        	  }
          tLock.unlock();
        }
    	printf("made it out of eventhappening loop. index is %d\n", index);

        }  // for loop

    	fclose(fp);  // timing.txt file.


        // trigger any remaining events in Kassiopeia so that its thread can finish.
    	while (fRunInProgress)
    		if ((fRunInProgress) && (fKassEventReady)) WakeBeforeEvent();





//    printf("for loop is ended\n");

        Kassiopeia.join();  // wait for Kassiopeia to finish.

//   printf("kass thread has joined.\n");

        return true;
    }

} /* namespace locust */
