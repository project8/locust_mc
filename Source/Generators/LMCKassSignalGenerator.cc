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
static std::string gxml_filename = "blank.xml";

FILE *fp2 = fopen("modeexctiation.txt","wb");  // time stamp checking.
FILE *fp3 = fopen("fabsfakemodeexctiation.txt","wb");  // time stamp checking.


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

    bool KassSignalGenerator::Configure( const scarab::param_node* aParam )
    {
        if( aParam == NULL) return true;
        fLO_Frequency = LO_FREQUENCY;

        if( aParam->has( "lo-frequency" ) )
        {
            fLO_Frequency = aParam->get_value< double >( "lo-frequency" );
        }

        if( aParam->has( "xml-filename" ) )
        {
            gxml_filename = aParam->get_value< std::string >( "xml-filename" );
        }

        return true;
    }

    void KassSignalGenerator::Accept( GeneratorVisitor* aVisitor ) const
    {
        aVisitor->Visit( this );
        return;
    }


  static void* KassiopeiaInit()
    {
        //cout << gxml_filename; getchar();
        const std::string & afile = gxml_filename;
    	RunKassiopeia *RunKassiopeia1 = new RunKassiopeia;
    	RunKassiopeia1->Run(afile);
    	delete RunKassiopeia1;
    }


    static void WakeBeforeEvent()
    {

        fPreEventCondition.notify_one();
        return;
    }


    static bool ReceivedKassReady()
    {

        if( !fKassEventReady)
        {
            std::unique_lock< std::mutex >tLock( fKassReadyMutex );
            fKassReadyCondition.wait( tLock );
            printf("LMC Got the fKassReadyCondition signal\n");
        }

        if (fFalseStartKassiopeia)  // workaround for some Macs
        {
            std::unique_lock< std::mutex >tLock( fKassReadyMutex );
            fKassReadyCondition.wait( tLock );
        }

        return true;
    }


    void* KassSignalGenerator::FilterNegativeFrequenciesNew(Signal* aSignal) const
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
            // Construct complex voltage.
            for( unsigned index = 0; index < windowsize; ++index )
            {
                SignalComplex[index][0] = aSignal->LongSignalTimeNew()[ nwin*windowsize + index ][0];
                SignalComplex[index][1] = aSignal->LongSignalTimeNew()[ nwin*windowsize + index ][1];
                //if (index==20000) {printf("signal 20000 is %g\n", aSignal->SignalTime()[index]); getchar();}
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
                // normalize and take the real part of the reverse transform, for digitization.
                //aSignal->SignalTime()[ nwin*windowsize + index ] = SignalComplex[index][0]/norm;
                aSignal->LongSignalTimeNew()[ nwin*windowsize + index ][0] = SignalComplex[index][0]/norm;
                //if (index>=20000) {printf("filtered signal is %g\n", aSignal->SignalTime()[index]); getchar();}
            }

        }




        delete SignalComplex;
        delete FFTComplex;

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
            // Construct complex voltage.
            for( unsigned index = 0; index < windowsize; ++index )
            {
                SignalComplex[index][0] = aLongSignal[ nwin*windowsize + index ];
                SignalComplex[index][1] = ImaginarySignal[ nwin*windowsize + index ];
                //if (index==20000) {printf("signal 20000 is %g\n", aSignal->SignalTime()[index]); getchar();}
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
                // normalize and take the real part of the reverse transform, for digitization.
                //aSignal->SignalTime()[ nwin*windowsize + index ] = SignalComplex[index][0]/norm;
                aLongSignal[ nwin*windowsize + index ] = SignalComplex[index][0]/norm;
                //if (index>=20000) {printf("filtered signal is %g\n", aSignal->SignalTime()[index]); getchar();}
            }

        }


        delete SignalComplex;
        delete FFTComplex;

    }

    void* KassSignalGenerator::DriveAntenna(int PreEventCounter, unsigned index, Signal* aSignal, double* anImaginarySignal) const
    {
        double tDopplerFrequencyAntenna = 0.;  // Doppler shifted cyclotron frequency in Hz.
        double tDopplerFrequencyShort = 0.;  
        double RealVoltage1 = 0.;
        double ImagVoltage1 = 0.;
        double RealVoltage2 = 0.;
        double ImagVoltage2 = 0.;
        double tCutOffFrequency = 2. * KConst::Pi() * KConst::C() * 1.841 / 2. / PI / 0.00502920; // rad/s, TE11

        locust::Particle tParticle = fParticleHistory.back();

        //Set as positive, even though really negative
        double tLarmorPower = tParticle.GetLarmorPower();
        double tCyclotronFrequency = tParticle.GetCyclotronFrequency()/2./KConst::Pi();
        if (tCyclotronFrequency < 15.e9) {printf("check 2PI in fcyc\n"); getchar();}
        double tVelocityZ = tParticle.GetVelocity().Z();
        double tGroupVelocity = KConst::C() * sqrt( 1. - pow(tCutOffFrequency/( 2.*KConst::Pi()*tCyclotronFrequency  ), 2.) );
        double tGammaZ = 1. / sqrt( 1. - pow(tVelocityZ / tGroupVelocity , 2. ) ); //generalization of lorentz factor to XXX mode waveguides, using only axial velocity of electrons

        //printf("paused in Locust! zvelocity is %g\n", zvelocity); getchar();

//        printf("GroupVelocity is %g, tCutOffFreq is %g, 2PIfcyc is %g\n", tGroupVelocity, tCutOffFrequency, 2.*KConst::Pi()*tCyclotronFrequency); getchar();
        tDopplerFrequencyAntenna = tCyclotronFrequency * tGammaZ *( 1. - tVelocityZ / tGroupVelocity);
        tDopplerFrequencyShort = tCyclotronFrequency *  tGammaZ *( 1. + tVelocityZ / tGroupVelocity);

        double tPositionZ = tParticle.GetPosition().Z();

        if (PreEventCounter > 0)
        {
            // initialize phases.
            phi_t1 = 2.*KConst::Pi()*(CENTER_TO_ANTENNA - tPositionZ) / (tGroupVelocity / tDopplerFrequencyAntenna);
            phi_t2 = 2.*KConst::Pi()*(tPositionZ + 2.*CENTER_TO_SHORT + CENTER_TO_ANTENNA) / (tGroupVelocity / tDopplerFrequencyShort);
            phi_shortTM01[0] = 0.;  // this gets advanced in the step modifier.
            phi_polarizerTM01[0] = 0.;  // this gets advanced in the step modifier.
        }

        //printf("PreEventCounter is %d and phi_t1 is %f and phi_t2 is %f\n", PreEventCounter, phi_t1, phi_t2); getchar();

        phi_t1 += 2.*KConst::Pi()*tDopplerFrequencyAntenna * fDigitizerTimeStep;
        phi_t2 += 2.*KConst::Pi()*tDopplerFrequencyShort * fDigitizerTimeStep;
        phiLO_t += 2.* KConst::Pi() * fLO_Frequency * fDigitizerTimeStep;
        RealVoltage1 = cos( phi_t1 - phiLO_t ); // + cos( phi_t1 + phiLO_t ));
        ImagVoltage1 = sin( phi_t1 - phiLO_t ); // + cos( phi_t1 + phiLO_t - PI/2.));
        RealVoltage2 = cos( phi_t2 - phiLO_t ); // + cos( phi_t2 + phiLO_t ));
        ImagVoltage2 = sin( phi_t2 - phiLO_t ); // + cos( phi_t2 + phiLO_t - PI/2.));

        //RealVoltage2 = 0.;  // take out short.
        //ImagVoltage2 = 0.;  // take out short.
        
        aLongSignal[ index ] += TE11ModeExcitation() * sqrt(tLarmorPower) * (RealVoltage1 + RealVoltage2);
        anImaginarySignal[ index ] += TE11ModeExcitation() * sqrt(tLarmorPower) * (ImagVoltage1 + ImagVoltage2);
        aSignal->LongSignalTimeNew()[ index ][0] += TE11ModeExcitation() * sqrt(tLarmorPower) * (RealVoltage1 + RealVoltage2);
        aSignal->LongSignalTimeNew()[ index ][1] += TE11ModeExcitation() * sqrt(tLarmorPower) * (ImagVoltage1 + ImagVoltage2);

        //if (t_old > 0.004)
/*
        {
            printf("driving antenna, ModeExcitation is %g\n\n", TE11ModeExcitation());
            printf("Realvoltage1 is %g and Realvoltage2 is %g\n", RealVoltage1, RealVoltage2);
            printf("Locust says:  signal %d is %g and zposition is %g and zvelocity is %g and sqrtLarmorPower is %g and "
            		"  fcyc is %.10g and tDopplerFrequency is %g and GammaZ is %.10g\n\n\n",
            index, aLongSignal[ index ], tPositionZ, tVelocityZ, pow(tLarmorPower,0.5), tCyclotronFrequency, tDopplerFrequencyAntenna, tGammaZ);
            getchar();
        }

        printf("fLO_Frequency is %g\n", fLO_Frequency); getchar();
*/

        t_old += fDigitizerTimeStep;  // advance time here instead of in step modifier.  This preserves the freefield sampling.

	  
    }

    double KassSignalGenerator::TE11ModeExcitation() const
    {
        double kc = 1.841/0.00502920;
        locust::Particle tParticle = fParticleHistory.back();
        double tPositionX = tParticle.GetPosition().X();
        double tPositionY = tParticle.GetPosition().Y();
        double r = sqrt( tPositionX*tPositionX + tPositionY*tPositionY);

        // fraction of emitted power that goes into TE11.
        // XXX
        double tCoupling = 119116./168.2 * 2./KConst::Pi() * 4./(2. * KConst::Pi()) / kc/2. * ( (j0(kc*r) - jn(2,kc*r)) +
                (j0(kc*r) + jn(2, kc*r)) );

        return sqrt(tCoupling);  // field amplitude is sqrt of power going into field.
    }


    bool KassSignalGenerator::DoGenerate( Signal* aSignal ) const
    {
        // temporary IQ patch.  Define and initialize ImaginarySignal.
        double *ImaginarySignal = new double[10*aSignal->TimeSize()];
        for( unsigned index = 0; index < 10*aSignal->TimeSize(); ++index )
        {
            ImaginarySignal[ index ] = 0.;
            aLongSignal[ index ] = 0.;  // long record for oversampling.
        }

        //n samples for event spacing.
        int PreEventCounter = 0;
        int NPreEventSamples = 150000;
        fPhaseIISimulation = true;

        //FILE *fp = fopen("timing.txt","wb");  // time stamp checking.
        //fprintf(fp, "testing\n");

        std::thread Kassiopeia (KassiopeiaInit);     // spawn new thread
        fRunInProgress = true;
        fKassEventReady = false;

        for( unsigned index = 0; index < 10*aSignal->TimeSize(); ++index )
        {
            if ((!fEventInProgress) && (fRunInProgress) && (!fPreEventInProgress))
            {
                if (ReceivedKassReady()) fPreEventInProgress = true;
                printf("LMC says it ReceivedKassReady()\n");
            }

            if (fPreEventInProgress)
            {
                PreEventCounter += 1;
                //printf("preeventcounter is %d\n", PreEventCounter);

                if (PreEventCounter > NPreEventSamples)  // finished noise samples.  Start event.
                {
                    fPreEventInProgress = false;  // reset.
                    fEventInProgress = true;
                    printf("LMC about to WakeBeforeEvent()\n");
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

        FilterNegativeFrequenciesNew(aSignal);
//        FilterNegativeFrequencies(aSignal, ImaginarySignal);
//        delete ImaginarySignal;

        //fclose(fp);  // timing.txt file.
        //fclose(fp2); // mode excitation file.

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
