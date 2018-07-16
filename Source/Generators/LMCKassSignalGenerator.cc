/*
 * LMCKassSignalGenerator.cc
 *
 *  Created on: Mar 12, 2014
 *      Author: nsoblath
 */

#include "LMCKassSignalGenerator.hh"
#include "LMCEventHold.hh"
#include "LMCRunKassiopeia.hh"
#include "LMCSimulationController.hh"
#include <chrono>


#include "logger.hh"
#include <thread>

#include "LMCGlobalsDeclaration.hh"

namespace locust
{
    LOGGER( lmclog, "KassSignalGenerator" );

    MT_REGISTER_GENERATOR(KassSignalGenerator, "kass-signal");

    KassSignalGenerator::KassSignalGenerator( const std::string& aName ) :
            Generator( aName ),
            fLO_Frequency( 0.),
            gxml_filename("blank.xml"),
            gpitchangle_filename("blank.xml"),
            phi_t1(0.),
            phi_t2(0.),
            phiLO_t(0.),
            EventStartTime(-99.),
            EventToFile(0)

    {
        fRequiredSignalState = Signal::kTime;
    }

    KassSignalGenerator::~KassSignalGenerator()
    {
    }

    bool KassSignalGenerator::Configure( const scarab::param_node* aParam )
    {
        if( aParam == NULL) return true;

        if( aParam->has( "lo-frequency" ) )
        {
            fLO_Frequency = aParam->get_value< double >( "lo-frequency" );
        }

        if( aParam->has( "xml-filename" ) )
        {
            gxml_filename = aParam->get_value< std::string >( "xml-filename" );
            gpitchangle_filename = aParam->get_value< std::string >( "pitchangle-filename" );
        }



        return true;
    }

    void KassSignalGenerator::Accept( GeneratorVisitor* aVisitor ) const
    {
        aVisitor->Visit( this );
        return;
    }


  static void* KassiopeiaInit(const std::string &aFile)
    {
    	RunKassiopeia *RunKassiopeia1 = new RunKassiopeia;
    	RunKassiopeia1->Run(aFile);
    	delete RunKassiopeia1;

        return 0;
    }


    static void WakeBeforeEvent()
    {

        fPreEventCondition.notify_one();
        return;
    }


    static bool ReceivedKassReady()
    {
    	std::this_thread::sleep_for(std::chrono::milliseconds(100));
		printf("LMC about to wait ..\n");

        if( !fKassEventReady)
        {
            std::unique_lock< std::mutex >tLock( fKassReadyMutex );
            fKassReadyCondition.wait( tLock );
        }

        if (fFalseStartKassiopeia)  // workaround for some Macs
        {
            std::unique_lock< std::mutex >tLock( fKassReadyMutex );
            fKassReadyCondition.wait( tLock );

        }

        return true;
    }

    void* KassSignalGenerator::DriveAntenna(int PreEventCounter, unsigned index, Signal* aSignal, FILE *fp)
    {
      double gain=1.;  //voltage gain
        double tDopplerFrequencyAntenna = 0.;  // Doppler shifted cyclotron frequency in Hz.
        double tDopplerFrequencyShort = 0.;  
        double RealVoltage1 = 0.;
        double ImagVoltage1 = 0.;
        double RealVoltage2 = 0.;
        double ImagVoltage2 = 0.;
        double tCutOffFrequency = 0.;
        if (Project8Phase == 2)
          {
          tCutOffFrequency = 2. * LMCConst::C() * 1.841 /( 2. * 0.00502920); // rad/s, TE11
          }
        else if (Project8Phase == 1)
          {
          tCutOffFrequency = LMCConst::C() * LMCConst::Pi() / 10.668e-3; // a in m
          }

        locust::Particle tParticle = fParticleHistory.back();
        RunLengthCalculator RunLengthCalculator1;


        //Set as positive, even though really negative for the particle.
        double tLarmorPower = tParticle.GetLarmorPower();
        double tCyclotronFrequency = tParticle.GetCyclotronFrequency()/2./LMCConst::Pi();
        if (tCyclotronFrequency < 15.e9) {printf("check 2PI in fcyc\n"); getchar();}
        double tVelocityZ = tParticle.GetVelocity().Z();
        double tPitchAngle = tParticle.GetPitchAngle();
        double tGroupVelocity = LMCConst::C() * sqrt( 1. - pow(tCutOffFrequency/( 2.*LMCConst::Pi()*tCyclotronFrequency  ), 2.) );
        double tGammaZ = 1. / sqrt( 1. - pow(tVelocityZ / tGroupVelocity , 2. ) ); //generalization of lorentz factor to XXX mode waveguides, using only axial velocity of electrons

        //printf("paused in Locust! zvelocity is %g\n", zvelocity); getchar();

//        printf("GroupVelocity is %g, tCutOffFreq is %g, 2PIfcyc is %g\n", tGroupVelocity, tCutOffFrequency, 2.*LMCConst::Pi()*tCyclotronFrequency); getchar();
        tDopplerFrequencyAntenna = tCyclotronFrequency * tGammaZ *( 1. - tVelocityZ / tGroupVelocity);
        tDopplerFrequencyShort = tCyclotronFrequency *  tGammaZ *( 1. + tVelocityZ / tGroupVelocity);

        double tPositionZ = tParticle.GetPosition().Z();

        if (PreEventCounter > 0)
        {
            // initialize phases.
            phi_t1 = 2.*LMCConst::Pi()*(CENTER_TO_ANTENNA - tPositionZ) / (tGroupVelocity / tDopplerFrequencyAntenna);
	    //            printf("center_to_antenna is %f and tPositionZ is %f\n", CENTER_TO_ANTENNA, tPositionZ);
            if (Project8Phase==1)
              {
              phi_t2 = LMCConst::Pi()/2. + 2.*LMCConst::Pi()*(CENTER_TO_SHORT + CENTER_TO_ANTENNA) / 
                (tGroupVelocity / tDopplerFrequencyShort);  // phase of reflected field at antenna.
              }
            if (Project8Phase==2)
	      {
	      phi_t2 = 2.*LMCConst::Pi()*(tPositionZ + 2.*CENTER_TO_SHORT + CENTER_TO_ANTENNA) / 
                (tGroupVelocity / tDopplerFrequencyShort);  // terminator in
	      }
            EventStartTime = (double)index/RunLengthCalculator1.GetAcquisitionRate()/1.e6/aSignal->DecimationFactor();
            EventToFile = false;
        }


        if ((tPitchAngle>0.)&&(EventToFile==false))
          {
	    fprintf(fp, "%10.4g   %g\n", EventStartTime, tPitchAngle);
          EventToFile = true;
          }



        phi_t1 += 2.*LMCConst::Pi()*tDopplerFrequencyAntenna * fDigitizerTimeStep;
        phi_t2 += 2.*LMCConst::Pi()*tDopplerFrequencyShort * fDigitizerTimeStep;
        phiLO_t += 2.* LMCConst::Pi() * fLO_Frequency * fDigitizerTimeStep;
        RealVoltage1 = cos( phi_t1 - phiLO_t ); // + cos( phi_t1 + phiLO_t ));
        ImagVoltage1 = sin( phi_t1 - phiLO_t ); // + cos( phi_t1 + phiLO_t - PI/2.));
        RealVoltage2 = cos( phi_t2 - phiLO_t ); // + cos( phi_t2 + phiLO_t ));
        ImagVoltage2 = sin( phi_t2 - phiLO_t ); // + cos( phi_t2 + phiLO_t - PI/2.));

        
        if (Project8Phase == 2)
          {
	    RealVoltage2 *= 0.03;  // replace short with terminator.                        
	    ImagVoltage2 *= 0.03;  // replace short with terminator.                                          
          aSignal->LongSignalTimeComplex()[ index ][0] += gain*sqrt(50.)*TE11ModeExcitation() * sqrt(tLarmorPower) * (RealVoltage1 + RealVoltage2);
          aSignal->LongSignalTimeComplex()[ index ][1] += gain*sqrt(50.)*TE11ModeExcitation() * sqrt(tLarmorPower) * (ImagVoltage1 + ImagVoltage2);
          }
        else if (Project8Phase == 1)
          {  // assume 50 ohm impedance

	    //	    RealVoltage2 *= 0.25; // some loss at short.
	    aSignal->LongSignalTimeComplex()[ index ][0] += gain*sqrt(50.) * TE01ModeExcitation() * sqrt(tLarmorPower) * (RealVoltage1 + RealVoltage2);
	    aSignal->LongSignalTimeComplex()[ index ][1] += gain*sqrt(50.) * TE01ModeExcitation() * sqrt(tLarmorPower) * (ImagVoltage1 + ImagVoltage2);
	    

          }


	/*						
            printf("driving antenna, ModeExcitation is %g\n\n", TE11ModeExcitation());
            printf("Realvoltage1 is %g and Realvoltage2 is %g\n", RealVoltage1, RealVoltage2);
            printf("IMagVoltage1 is %g and ImagVoltage2 is %g\n", ImagVoltage1, ImagVoltage2);
            printf("Locust says:  signal %d is %g and zposition is %g and zvelocity is %g and sqrtLarmorPower is %g and "
            		"  fcyc is %.10g and tDopplerFrequency is %g and GammaZ is %.10g\n\n\n",
            index, aSignal->LongSignalTimeComplex()[ index ][0], tPositionZ, tVelocityZ, pow(tLarmorPower,0.5), tCyclotronFrequency, tDopplerFrequencyAntenna, tGammaZ);
//            getchar();


        printf("fLO_Frequency is %g\n", fLO_Frequency); getchar();
*/		
	

        t_old += fDigitizerTimeStep;  // advance time here instead of in step modifier.  This preserves the freefield sampling.
	  
        return 0;
    }

    double KassSignalGenerator::TE11ModeExcitation() const
    {
        double kc = 1.841/0.00502920;
        locust::Particle tParticle = fParticleHistory.back();
        double tPositionX = tParticle.GetPosition().X();
        double tPositionY = tParticle.GetPosition().Y();
        double r = sqrt( tPositionX*tPositionX + tPositionY*tPositionY);

        // field amplitude in TE11.  Apply this factor to Larmor Power and propagate it to the amplifier.
        double tCoupling = 119116./168.2 * 2./LMCConst::Pi() * 4./(2. * LMCConst::Pi()) / kc/2. * ( (j0(kc*r) - jn(2,kc*r)) +
                (j0(kc*r) + jn(2, kc*r)) );

        return tCoupling;  // field amplitude is sqrt of power going into field.
    }


    double KassSignalGenerator::TE01ModeExcitation() const
    {
   	    double dim1_wr42 = 10.668e-3; // a in m
        locust::Particle tParticle = fParticleHistory.back();
	    double tPositionX = tParticle.GetPosition().X() + dim1_wr42/2.;
	    double coupling = 0.63*sin(LMCConst::Pi()*tPositionX/dim1_wr42);  // avg over cyclotron orbit.
	    return coupling;
    }



    bool KassSignalGenerator::DoGenerate( Signal* aSignal )
    {
        //n samples for event spacing.
        int PreEventCounter = 0;
        int NPreEventSamples = 1500000;

        //FILE *fp = fopen("timing.txt","wb");  // time stamp checking.
        //fprintf(fp, "testing\n");
        FILE *fp = fopen(gpitchangle_filename.c_str(), "w");


        std::thread Kassiopeia (KassiopeiaInit,gxml_filename);     // spawn new thread
        fRunInProgress = true;
        fKassEventReady = false;
        int StartEventTimer = 0;

        for( unsigned index = 0; index < aSignal->DecimationFactor()*aSignal->TimeSize(); ++index )
        {
            if ((!fEventInProgress) && (fRunInProgress) && (!fPreEventInProgress))
            {
                if (ReceivedKassReady()) fPreEventInProgress = true;
                printf("LMC says it ReceivedKassReady()\n");                
		if (Project8Phase == 2)
                  {
		    if ((index-StartEventTimer) > 1000) 
		    {
                    break;
		    }
                  }
            }

            if (fPreEventInProgress)
            {
	        PreEventCounter += 1;

                if (PreEventCounter > NPreEventSamples)  // finished noise samples.  Start event.
                {
                    fPreEventInProgress = false;  // reset.
                    fEventInProgress = true;
                    printf("LMC about to WakeBeforeEvent()\n");
                    StartEventTimer = index;
                    WakeBeforeEvent();  // trigger Kass event.
                }
            }

            if (fEventInProgress)  // fEventInProgress
            if (fEventInProgress)  // check again.
            {
                std::unique_lock< std::mutex >tLock( fMutexDigitizer, std::defer_lock );
                tLock.lock();
                fDigitizerCondition.wait( tLock );
                if (fEventInProgress)
                {
                    DriveAntenna(PreEventCounter, index, aSignal, fp);
                    PreEventCounter = 0; // reset
                }
                tLock.unlock();
            }
        }  // for loop

        printf("finished signal loop.\n");


        fclose(fp);
	fRunInProgress = false;  // tell Kassiopeia to finish.
        fDoneWithSignalGeneration = true;  // tell LMCCyclotronRadExtractor
	//	if (fEventInProgress)
	//  if (ReceivedKassReady())
	    WakeBeforeEvent();
			     
         
        Kassiopeia.join();  // wait for Kassiopeia to finish.
	

        return true;
    }

} /* namespace locust */
