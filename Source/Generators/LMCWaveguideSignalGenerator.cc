/*
 * LMCWaveguideSignalGenerator.cc
 *
 *  Created on: Mar 12, 2014
 *      Author: nsoblath
 */

#include "LMCWaveguideSignalGenerator.hh"

#include "LMCRunKassiopeia.hh"

#include "logger.hh"

#include <chrono>
#include <thread>

namespace locust
{
    LOGGER( lmclog, "WaveguideSignalGenerator" );

    MT_REGISTER_GENERATOR(WaveguideSignalGenerator, "waveguide-signal");

    WaveguideSignalGenerator::WaveguideSignalGenerator( const std::string& aName ) :
    		Generator( aName ),
			fLO_Frequency( 0.),
			gxml_filename("blank.xml"),
			gpitchangle_filename("blank.xml"),
			fTruth( 0 ),
			fPhi_t1(0.),
			fPhi_t2(0.),
			fPhiLO_t(0.),
			fNPreEventSamples( 15 ),
			fEventStartTime(-99.),
			fEventToFile(0),
			fInterface( new KassLocustInterface() )

    {
        fRequiredSignalState = Signal::kTime;

        KLInterfaceBootstrapper::get_instance()->SetInterface( fInterface );
    }

    WaveguideSignalGenerator::~WaveguideSignalGenerator()
    {
    }

    bool WaveguideSignalGenerator::Configure( const scarab::param_node& aParam )
    {
        if( aParam.has( "lo-frequency" ) )
        {
            fLO_Frequency = aParam["lo-frequency"]().as_double();
        }

        if( aParam.has( "xml-filename" ) )
        {
            gxml_filename = aParam["xml-filename"]().as_string();
        }

        if( aParam.has( "pitchangle-filename" ) )
        {
            gpitchangle_filename = aParam["pitchangle-filename"]().as_string();
        }

        if( aParam.has( "truth" ) )
        {
            fTruth = aParam["truth"]().as_bool();
        }

        return true;
    }

    void WaveguideSignalGenerator::Accept( GeneratorVisitor* aVisitor ) const
    {
        aVisitor->Visit( this );
        return;
    }



    void WaveguideSignalGenerator::KassiopeiaInit(const std::string &aFile)
    {
        RunKassiopeia tRunKassiopeia;
        tRunKassiopeia.Run(aFile, fInterface);
        return;
    }


    void WaveguideSignalGenerator::WakeBeforeEvent()
    {
    	fInterface->fRunInProgress = false;
        fInterface->fPreEventCondition.notify_one();
        return;
    }


    bool WaveguideSignalGenerator::ReceivedKassReady()
    {
        std::this_thread::sleep_for(std::chrono::milliseconds(100));

        printf("LMC about to wait fKassEventReady is %d, fRunInProgress is %d\n ", fInterface->fKassEventReady, fInterface->fRunInProgress);

        if((fInterface->fRunInProgress)&&(!fInterface->fKassEventReady))
        {
            std::unique_lock< std::mutex >tLock( fInterface->fKassReadyMutex );
            fInterface->fKassReadyCondition.wait( tLock );
            return true;
        }
        else if (!fInterface->fRunInProgress)
        {
        	return false;
        }
        else if (fInterface->fKassEventReady)
        {
        	return true;
        }
        else
        {
            printf("I am stuck.\n"); getchar();
            return true;
        }


    }

    void* WaveguideSignalGenerator::DriveAntenna(int PreEventCounter, unsigned index, Signal* aSignal, FILE *fp)
    {

        double tDopplerFrequencyAntenna = 0.;  // Doppler shifted cyclotron frequency in Hz.
        double tDopplerFrequencyShort = 0.;  
        double RealVoltage1 = 0.;
        double ImagVoltage1 = 0.;
        double RealVoltage2 = 0.;
        double ImagVoltage2 = 0.;
        double tCutOffFrequency = 0.;  // Hz, placeholder

        tCutOffFrequency = LMCConst::C() * LMCConst::Pi() / 10.668e-3; // a in m


        int currentIndex = FindNode(fInterface->fTOld);
        locust::Particle tParticle = fInterface->fParticleHistory[currentIndex];
        tParticle.Interpolate(fInterface->fTOld);

        //Set as positive, even though really negative for the particle.
        double tLarmorPower = tParticle.GetLarmorPower();
        double tCyclotronFrequency = tParticle.GetCyclotronFrequency()/2./LMCConst::Pi();
        if (tCyclotronFrequency < 15.e9) {printf("check 2PI in fcyc\n"); getchar();}
        double tVelocityZ = tParticle.GetVelocity().Z();
        double tPitchAngle = tParticle.GetPitchAngle();
        double tGroupVelocity = LMCConst::C() * sqrt( 1. - pow(tCutOffFrequency/( 2.*LMCConst::Pi()*tCyclotronFrequency  ), 2.) );
        double tGammaZ = 1. / sqrt( 1. - pow(tVelocityZ / tGroupVelocity , 2. ) ); //generalization of lorentz factor to XXX mode waveguides, using only axial velocity of electrons

        tDopplerFrequencyAntenna = tCyclotronFrequency * tGammaZ *( 1. - tVelocityZ / tGroupVelocity);
        tDopplerFrequencyShort = tCyclotronFrequency *  tGammaZ *( 1. + tVelocityZ / tGroupVelocity);

        double tPositionZ = tParticle.GetPosition().Z();

        if (PreEventCounter > 0)
        {
            // initialize radiation phases.

            fPhi_t1 = 2.*LMCConst::Pi()*(fInterface->fCENTER_TO_ANTENNA - tPositionZ) / (tGroupVelocity / tDopplerFrequencyAntenna);

            fPhi_t2 = LMCConst::Pi()/2. + 2.*LMCConst::Pi()*(fInterface->fCENTER_TO_SHORT + fInterface->fCENTER_TO_ANTENNA) /
                    (tGroupVelocity / tDopplerFrequencyShort);  // phase of reflected field at antenna.

            fEventStartTime = (double)index/fAcquisitionRate/1.e6/aSignal->DecimationFactor();

            fEventToFile = false;
        }


        if ((tPitchAngle>0.)&&(fEventToFile==false))
        {
            fprintf(fp, "%10.4g   %g\n", fEventStartTime, tPitchAngle);
            fEventToFile = true;
        }


        fPhi_t1 += 2.*LMCConst::Pi()*tDopplerFrequencyAntenna * 1./(fAcquisitionRate*1.e6*aSignal->DecimationFactor());
        fPhi_t2 += 2.*LMCConst::Pi()*tDopplerFrequencyShort * 1./(fAcquisitionRate*1.e6*aSignal->DecimationFactor());
        fPhiLO_t += 2.* LMCConst::Pi() * fLO_Frequency * 1./(fAcquisitionRate*1.e6*aSignal->DecimationFactor());
        RealVoltage1 = cos( fPhi_t1 - fPhiLO_t ); // + cos( phi_t1 + phiLO_t ));  // antenna
        ImagVoltage1 = sin( fPhi_t1 - fPhiLO_t ); // + cos( phi_t1 + phiLO_t - PI/2.));
        RealVoltage2 = cos( fPhi_t2 - fPhiLO_t ); // + cos( phi_t2 + phiLO_t ));  // short
        ImagVoltage2 = sin( fPhi_t2 - fPhiLO_t ); // + cos( phi_t2 + phiLO_t - PI/2.));


        RealVoltage2 *= 0.03;  // replace short with terminator.
        ImagVoltage2 *= 0.03;  // replace short with terminator.
        aSignal->LongSignalTimeComplex()[ index ][0] += sqrt(50.)*ModeExcitation() * sqrt(tLarmorPower/2.) * (RealVoltage1 + RealVoltage2);
        aSignal->LongSignalTimeComplex()[ index ][1] += sqrt(50.)*ModeExcitation() * sqrt(tLarmorPower/2.) * (ImagVoltage1 + ImagVoltage2);


        fInterface->fTOld += 1./(fAcquisitionRate*1.e6*aSignal->DecimationFactor());  // advance time here instead of in step modifier.  This preserves the freefield sampling.

        return 0;
    }

    //Return index of fParticleHistory particle closest to the time we are evaluating
    int WaveguideSignalGenerator::FindNode(double tNew) const
    {
        std::deque<locust::Particle>::iterator it;
        it = std::upper_bound( fInterface->fParticleHistory.begin() , fInterface->fParticleHistory.end() , tNew, [] (const double &a , const locust::Particle &b) { return a < b.GetTime();} );

        int tNodeIndex = it - fInterface->fParticleHistory.begin();

        return tNodeIndex;
    }


    double WaveguideSignalGenerator::ModeExcitation() const
    {
    	locust::Particle tParticle = fInterface->fParticleHistory.back();
        double tPositionX = tParticle.GetPosition().X();
        double tPositionY = tParticle.GetPosition().Y();
        double r = sqrt( tPositionX*tPositionX + tPositionY*tPositionY);
        double tCoupling = 1.0; // placeholder
        return tCoupling;  // field amplitude is sqrt of power fraction.
    }


    bool WaveguideSignalGenerator::DoGenerate( Signal* aSignal )
    {
        int PreEventCounter = 0;

        FILE *fp = fopen(gpitchangle_filename.c_str(), "w");

        fInterface->fKassTimeStep = 1./(fAcquisitionRate*1.e6*aSignal->DecimationFactor());
        std::thread tKassiopeia (&WaveguideSignalGenerator::KassiopeiaInit, this, gxml_filename);     // spawn new thread

        for( unsigned index = 0; index < aSignal->DecimationFactor()*aSignal->TimeSize(); ++index )
        {
            if ((!fInterface->fEventInProgress) /*&& (fInterface->fRunInProgress)*/ && (!fInterface->fPreEventInProgress))
            {
            	if (ReceivedKassReady()) fInterface->fPreEventInProgress = true;
            	else
            	{
            		printf("breaking\n");
            		break;
            	}

            	printf("LMC says it ReceivedKassReady(), fRunInProgress is %d\n", fInterface->fRunInProgress);

            }

            if ((fInterface->fPreEventInProgress)&&(fInterface->fRunInProgress))
            {
                PreEventCounter += 1;

                if (((!fTruth)&&(PreEventCounter > fNPreEventSamples))||((fTruth)&&(PreEventCounter > fNPreEventSamples)&&(index%(8192*aSignal->DecimationFactor())==0)  ))// finished pre-samples.  Start event.
                {
                    fInterface->fPreEventInProgress = false;  // reset.
                    fInterface->fEventInProgress = true;
                    printf("LMC about to WakeBeforeEvent()\n");
                    WakeBeforeEvent();  // trigger Kass event.
                }
            }

            if ((fInterface->fEventInProgress)&&(!fInterface->fKassEventReady))  // fEventInProgress
            {
                    std::unique_lock< std::mutex >tLock( fInterface->fMutexDigitizer, std::defer_lock );
                    tLock.lock();
                    fInterface->fDigitizerCondition.wait( tLock );
                    if (fInterface->fEventInProgress)
                    {
                        DriveAntenna(PreEventCounter, index, aSignal, fp);
                        PreEventCounter = 0; // reset
                    }
                    tLock.unlock();
            }
        }  // for loop


        fInterface->fDoneWithSignalGeneration = true;
        fclose(fp);
        printf("finished signal loop.\n");
        WakeBeforeEvent();
        tKassiopeia.join();  // finish thread

        return true;
    }

} /* namespace locust */
