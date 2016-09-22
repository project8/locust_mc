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
double dummy_time = 0.;  // time for checking mode excitation.

    	FILE *fp2 = fopen("modeexctiation.txt","wb");  // time stamp checking.
    	FILE *fp3 = fopen("fabsfakemodeexctiation.txt","wb");  // time stamp checking.


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
    	const string & afile = "/home/slocum/project8/Project8_by_Devin.xml";
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
    if( !fKassEventReady)  // check again.
    {
    std::unique_lock< std::mutex >tLock( fMutex );
    fKassReadyCondition.wait( tLock );
    printf("LMC Got the fKassReadyCondition signal\n");
    return true;
    }
    else {return true;}
}




void* KassSignalGenerator::FilterNegativeFrequencies(Signal* aSignal, double *ImaginarySignal) const
{

  int nwindows = 8;
  int windowsize = aSignal->TimeSize()/nwindows;


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


// Construct complex voltage.
    for( unsigned index = 0; index < windowsize; ++index )
    {
        SignalComplex[index][0] = aSignal->SignalTime()[ nwin*windowsize + index ];
        SignalComplex[index][1] = ImaginarySignal[ nwin*windowsize + index ];
//        if (index==20000) {printf("signal 20000 is %g\n", aSignal->SignalTime()[index]); getchar();}

    }

    fftw_execute(ForwardPlan);


    

//  Complex filter to set power at negative frequencies to 0.

    for( unsigned index = windowsize/2; index < windowsize; ++index )
    {
        FFTComplex[index][0] = 0.;
        FFTComplex[index][1] = 0.;
    }

    


    /*
    // Complex low pass filter 90% of Nyquist.
        for( unsigned index = 0; index < aSignal->TimeSize(); ++index )
        {
        	if ((index > aSignal->TimeSize()/2*90/100 ) && (index < aSignal->TimeSize()/2*110/100 ))
            FFTComplex[index][0] = 0.;
            FFTComplex[index][1] = 0.;
        }

    */    



    fftw_execute(ReversePlan);

    double norm = (double)(windowsize);

    for( unsigned index = 0; index < windowsize; ++index )
    {
    	// normalize and take the real part of the reverse transform, for digitization.
      aSignal->SignalTime()[ nwin*windowsize + index ] = SignalComplex[index][0]/norm;
      //    if (index>=20000) {printf("filtered signal is %g\n", aSignal->SignalTime()[index]); getchar();}
    }



    }


delete SignalComplex;
delete FFTComplex;


}




void* KassSignalGenerator::DriveAntenna(unsigned index, Signal* aSignal, double* ImaginarySignal) const
{

	double dt = 5.e-9; // seconds, this might need to come from Kassiopeia and be exact.
    double fprime = 0.;  // Doppler shifted cyclotron frequency in Hz.
    double RealVoltagePhase = 0.;
    double ImagVoltagePhase = 0.;
    double GroupVelocity = 0.;
    double SpeedOfLight = 2.99792458e8; // m/s
    double CutOffFrequency = SpeedOfLight * PI / 10.668e-3; // a in m         
                                

//                  printf("paused in Locust! zvelocity is %g\n", zvelocity); getchar();

    GroupVelocity = SpeedOfLight * pow( 1. - pow(CutOffFrequency/(2.*PI*fcyc), 2.) , 0.5);
    //    printf("GroupVelocity is %g, CutOffFreq is %g, 2PIfcyc is %g\n", GroupVelocity, CutOffFrequency, 2.*PI*fcyc); getchar();
           fprime = fcyc*GammaZ*(1.-zvelocity/GroupVelocity);
           phi_t += 2.*PI*fprime*dt;
           phiLO_t += 2.*PI*LO_FREQUENCY*dt;
           RealVoltagePhase = cos( phi_t - phiLO_t ); // + cos( phi_t + phiLO_t ));
           ImagVoltagePhase = cos( phi_t - phiLO_t - PI/2.); // + cos( phi_t + phiLO_t - PI/2.));


           aSignal->SignalTime()[ index ] += AverageModeExcitation()*pow(LarmorPower,0.5)*RealVoltagePhase;
           ImaginarySignal[ index ] += AverageModeExcitation()*pow(LarmorPower,0.5)*ImagVoltagePhase;


	   //	   	   	   	   	   	   	   	          printf("driving antenna, ModeExcitation is %g\n\n", AverageModeExcitation());
	   //	   	   	              printf("Locust says:  signal %d is %g and t is %g and zvelocity is %g and sqrtLarmorPower is %g and fcyc is %.10g and fprime is %g and GammaZ is %.10g\n",
	   //	   	                      index, aSignal->SignalTime()[ index ], t_poststep, zvelocity, pow(LarmorPower,0.5), fcyc, fprime, GammaZ);
	   //	   	                     getchar();


}



double KassSignalGenerator::ModeExcitation() const
{
    double dim1_wr42 = 10.668e-3; // a in m
    double dim2_wr42 = 4.318e-3; // b in m
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
//	printf("x is %f and Ey is %g and yvelocity is %g and xvelocity is %g and EdotV is %f\n", x, Ey, vy, vx, EdotV);

//	fprintf(fp2, "%.10g %.10g %.10g %.10g ", EdotV, vx, vy, X);  // checking mode excitation.


return EdotV;
}


double KassSignalGenerator::AverageModeExcitation() const
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

	// E dot v / (Emax v) * 0.637 / sqrt(2) for analytic avg coupling/orbit and half power lost in opposite direction.
	// or E dot v / (Emax v) * 0.637 for analytic avg coupling/orbit and reflecting short.
	double AverageEdotV = fabs(Ey/EyMax) *0.637;

//	printf("x is %f and Ey is %g and yvelocity is %g and xvelocity is %g and EdotV is %f\n", x, Ey, vy, vx, EdotV);

//	fprintf(fp2, "%.10g %.10g %.10g %.10g ", EdotV, vx, vy, X);  // checking mode excitation.

return AverageEdotV;
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



    bool KassSignalGenerator::DoGenerate( Signal* aSignal ) const
    {

      // temporary IQ patch.  Define and initialize ImaginarySignal.
    	double *ImaginarySignal = new double[aSignal->TimeSize()]; 
        for( unsigned index = 0; index < aSignal->TimeSize(); ++index )
          ImaginarySignal[ index ] = 0.;

//      n samples for event spacing.
    	int PreEventCounter = 0;
	int NPreEventSamples = 150000;

//    	FILE *fp = fopen("timing.txt","wb");  // time stamp checking.
//    	fprintf(fp, "testing\n");


    	std::thread Kassiopeia (KassiopeiaInit);     // spawn new thread
    	fRunInProgress = true;


	for( unsigned index = 0; index < aSignal->TimeSize(); ++index )


    	{

        if ((!fEventInProgress) && (fRunInProgress) && (!fPreEventInProgress))
        	{
        	if (ReceivedKassReady()) fPreEventInProgress = true;
		               
        	}


    	if (fPreEventInProgress)
    	  {
    	  PreEventCounter += 1;
	  //	            printf("preeventcounter is %d\n", PreEventCounter);

          if (PreEventCounter > NPreEventSamples)  // finished noise samples.  Start event.
        	  {
        	  fPreEventInProgress = false;  // reset.
                  fEventInProgress = true;
        	  PreEventCounter = 0;  // reset.
    		  phi_t = 0.; // initialize antenna voltage phase in radians.
  		  phiLO_t = 0.; // initialize voltage phase of LO;
                  printf("LMC about to wakebeforeevent\n");
	          WakeBeforeEvent();  // trigger Kass event.
        	  }
    	  }


    	if (fEventInProgress)  // fEventInProgress
	if (fEventInProgress)  // check again.
        {
//    	  printf("waiting for digitizer trigger ... index is %d\n", index);
          std::unique_lock< std::mutex >tLock( fMutexDigitizer, std::defer_lock );
          tLock.lock();
          fDigitizerCondition.wait( tLock );
          if (fEventInProgress)
        	  {
        	  DriveAntenna(index, aSignal, ImaginarySignal);
        	  }
          tLock.unlock();
        }

        }  // for loop


	FilterNegativeFrequencies(aSignal, ImaginarySignal);
        delete ImaginarySignal;


//    	fclose(fp);  // timing.txt file.
//        fclose(fp2); // mode excitation file.

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
