/*
 * LMCCyclotronRadiationExtractor.cc
 *
 *  Created on: Mar 13, 2016
 *      Author: nsoblath
 */

#include "LMCCyclotronRadiationExtractor.hh"
#include "LMCGlobalsDeclaration.hh"

using namespace Kassiopeia;
namespace locust
{

    CyclotronRadiationExtractor::CyclotronRadiationExtractor():
    		fP8Phase( 0 ),
    		fPitchAngle( -99. )
    {
    }

    CyclotronRadiationExtractor::CyclotronRadiationExtractor( const CyclotronRadiationExtractor& aOrig ):
    		fP8Phase( 0 ),
    		fPitchAngle( -99. )
    {
    }

    CyclotronRadiationExtractor::~CyclotronRadiationExtractor()
    {
    }

    CyclotronRadiationExtractor* CyclotronRadiationExtractor::Clone() const
    {
        return new CyclotronRadiationExtractor( *this );
    }


    void CyclotronRadiationExtractor::SetP8Phase (int P8Phase )
    {
    	fP8Phase = P8Phase;
    	Project8Phase = P8Phase;
        if (P8Phase==1)
          {
	    CENTER_TO_SHORT = 0.047; // m
            CENTER_TO_ANTENNA = 0.045; // m
          }
	if (P8Phase==2)
	  {
	    CENTER_TO_SHORT = 0.075; // m
	    CENTER_TO_ANTENNA = 0.075; // m  
          }
    }

    bool CyclotronRadiationExtractor::ExecutePreStepModification( KSParticle& anInitialParticle, KSParticleQueue& aQueue )
    {
    	return true;
    }

    double CyclotronRadiationExtractor::GetGroupVelocityTM01(KSParticle& aFinalParticle)
    {
        const double SpeedOfLight = LMCConst::C();
        double CutOffFrequency = 2. * LMCConst::Pi() * SpeedOfLight * 2.405 / 2. / LMCConst::Pi() / 0.00502920; // rad/s
        double cyclotronFrequency = aFinalParticle.GetCyclotronFrequency();
        double GroupVelocity = SpeedOfLight * sqrt( 1. - pow(CutOffFrequency/(2. * LMCConst::Pi() * cyclotronFrequency), 2.) );
    	return GroupVelocity;
    }

    double CyclotronRadiationExtractor::GetGroupVelocityTE11(KSParticle& aFinalParticle)
    {
        const double SpeedOfLight = LMCConst::C(); // m/s
        double CutOffFrequency = 2. * LMCConst::Pi() * SpeedOfLight * 1.841 / 2. / LMCConst::Pi() / 0.00502920; // rad/s
        double cyclotronFrequency = aFinalParticle.GetCyclotronFrequency();
        double GroupVelocity = SpeedOfLight * sqrt( 1. - pow(CutOffFrequency / (2. * LMCConst::Pi() * cyclotronFrequency), 2.) );
    	return GroupVelocity;
    }

    double CyclotronRadiationExtractor::GetGroupVelocityTE01(KSParticle& aFinalParticle)  // Phase 1
     {
       double SpeedOfLight = LMCConst::C(); // m/s
         double CutOffFrequency = SpeedOfLight * LMCConst::Pi() / 10.668e-3; // a in m
         double fcyc = aFinalParticle.GetCyclotronFrequency();
         double GroupVelocity = SpeedOfLight * pow( 1. - pow(CutOffFrequency/(2.*LMCConst::Pi()*fcyc), 2.) , 0.5);
 //        printf("GroupVelocity is %g\n", GroupVelocity); getchar();
     	return GroupVelocity;
     }



    double CyclotronRadiationExtractor::GetCouplingFactorTE11(KSParticle& aFinalParticle)
    {
    	double kc = 1.841/0.00502920;
    	double x = aFinalParticle.GetPosition().GetX();
    	double y = aFinalParticle.GetPosition().GetY();
        
    	double r = sqrt( x * x + y * y);
    	double coupling = 119116./168.2 * 2./LMCConst::Pi() * 4./(2.*LMCConst::Pi()) / kc/2. * ( (j0(kc*r) - jn(2,kc*r)) +
    			(j0(kc*r) + jn(2, kc*r)) );
    	return coupling*coupling;
    }

    double CyclotronRadiationExtractor::GetCouplingFactorTM01(KSParticle& aFinalParticle)
    {
    	double kc = 2.405/0.00502920;
    	double x = aFinalParticle.GetPosition().GetX();
    	double y = aFinalParticle.GetPosition().GetY();
    	double r = sqrt(x*x + y*y);
    	double coupling =   146876.5/168.2 * 2./LMCConst::Pi() * 4./(2.*LMCConst::Pi()) / kc * j1(kc*r);
    	return coupling*coupling;
    }

    double CyclotronRadiationExtractor::GetCouplingFactorTE01(KSParticle& aFinalParticle)  // Phase 1
    {
    	double dim1_wr42 = 10.668e-3; // a in m
    	double x = aFinalParticle.GetPosition().GetX() + dim1_wr42/2.;
    	double coupling = 0.63*sin(LMCConst::Pi()*x/dim1_wr42);  // avg over cyclotron orbit.
    	return coupling*coupling;
    }


    double CyclotronRadiationExtractor::GetTE01FieldAfterOneBounce(KSParticle& anInitialParticle, KSParticle& aFinalParticle)
    {
        double fcyc = aFinalParticle.GetCyclotronFrequency();
        double GroupVelocity = GetGroupVelocityTE01(aFinalParticle);
    	double zvelocity = aFinalParticle.GetVelocity().GetZ();
    	double zPosition = aFinalParticle.GetPosition().GetZ();
        double GammaZ = 1.0/pow(1.0-pow(zvelocity/GetGroupVelocityTE01(aFinalParticle),2.),0.5);

    	double fprime_short = fcyc*GammaZ*(1.+zvelocity/GroupVelocity);
    	double phi_shortTE01 = LMCConst::Pi()/2. + 2.*LMCConst::Pi()*(zPosition+CENTER_TO_SHORT)/(GroupVelocity/fprime_short);  // phase of reflected field at position of electron.
//        double FieldFromShort = cos(phi_shortTM01);  // no resonant enhancement.
        double FieldFromShort = cos(0.) + cos(phi_shortTE01); // yes resonant enhancement.
	//        printf("\n phi_shortTE01 is %f and phi_shortTE01/(2PI) is %f\n", phi_shortTE01, phi_shortTE01/2./LMCConst::Pi());
        //printf("phi_shortTE01 at antenna should be %f\n", phi_shortTE01+2.*LMCConst::Pi()*(0.045/(GroupVelocity/fprime_short)));
        //printf("phi_otherdirection should be %f\n", 2.*LMCConst::Pi()*(0.045/(GroupVelocity/fprime_short)));
	//printf("field after short is %f\n\n", FieldFromShort); getchar();
 

        return FieldFromShort;  // Phase 1

    }

    double CyclotronRadiationExtractor::GetTE11FieldAfterOneBounce(KSParticle& anInitialParticle, KSParticle& aFinalParticle)
    {
    	double dt = aFinalParticle.GetTime() - anInitialParticle.GetTime();
        double cyclotronFrequency = aFinalParticle.GetCyclotronFrequency();
        double GroupVelocity = GetGroupVelocityTE11(aFinalParticle);

    	double zVelocity = aFinalParticle.GetVelocity().GetZ();
    	double zPosition = aFinalParticle.GetPosition().GetZ();
        double GammaZ = 1.0/ pow( 1.0 - pow( zVelocity / GroupVelocity,2.),0.5);
    	double fprime_short = cyclotronFrequency*GammaZ*(1.+zVelocity / GroupVelocity);
    	double TE11FieldAfterOneBounce = 0.;
    	double phi_shortTE11 = 0.;

    	// Cu boundary condition gives PI/2 phase advancement to short.
        phi_shortTE11 = LMCConst::Pi()/2. + 2.*LMCConst::Pi()*(zPosition+CENTER_TO_SHORT)/(GroupVelocity/fprime_short);
        TE11FieldAfterOneBounce = cos(0.) + cos(phi_shortTE11);
        //printf("TE11FieldAfterOneBounce is %f\n", TE11FieldAfterOneBounce);

    	return TE11FieldAfterOneBounce;
    }

    double CyclotronRadiationExtractor::GetTM01FieldAfterBounces(KSParticle& anInitialParticle, KSParticle& aFinalParticle)
    {

    	double dt = aFinalParticle.GetTime() - anInitialParticle.GetTime();
        double tCyclotronFrequency = aFinalParticle.GetCyclotronFrequency();
        double GroupVelocity = GetGroupVelocityTM01(aFinalParticle);
        double tVelocityZ = aFinalParticle.GetVelocity().GetZ();
        double tPositionZ = aFinalParticle.GetPosition().GetZ();
        double GammaZ = 1.0/sqrt(1.0-pow(tVelocityZ / GetGroupVelocityTM01(aFinalParticle),2.) );

    	double fprime_short = tCyclotronFrequency*GammaZ*(1.+tVelocityZ/GroupVelocity);
//    	printf("tcyc is %g and GammaZ is %g and tVelocityZ is %g\n", tCyclotronFrequency, GammaZ, tVelocityZ);
    	double fprime_polarizer = tCyclotronFrequency*GammaZ*(1.-tVelocityZ/GroupVelocity);
    	double lambda_short = GroupVelocity/fprime_short;
//    	printf("GroupVelocity is %.10g and fprime_short is %.10g\n", GroupVelocity, fprime_short);
    	double lambda_polarizer = GroupVelocity/fprime_polarizer;
    	double FieldFromShort=0.;  // first doppler shift
    	double FieldFromPolarizer=0.; // other doppler shift
    	double TM01FieldAfterBounces = 0.;
    	int nbounces = 10;  // even number please.
        double phi_shortTM01 = 0.;
        double phi_polarizerTM01 = 0.;
        double dphi = 0.;

//   	       	printf("phi_shortTM01[0] is %.10g and LMCConst::Pi() is %.10g and z is %.10g and lambda is %.10g\n", phi_shortTM01[0], LMCConst::Pi(), tPositionZ, lambda_short);

            FieldFromShort = cos(0.); // starting field, after no reflections.  divide by 2 later.
            FieldFromPolarizer = cos(0.); // starting field, after no reflections.  divide by 2 later.


            for (int i=0; i<nbounces; i++)  // short-going wave, initially.  It already hit the short.
            {
            	if (i%2==0) // now toward polarizer.
            	  {
            	  phi_shortTM01 += 2.*LMCConst::Pi()*(CENTER_TO_SHORT + CENTER_TO_ANTENNA)/lambda_short + LMCConst::Pi();  // group velocity, with phase shift
            	  // BACK UP toward short, stop at electron and calculate field there:
            	  FieldFromShort += cos(phi_shortTM01 - 2.*LMCConst::Pi()*(CENTER_TO_ANTENNA - tPositionZ)/lambda_short);
            	  }
            	else  // now toward short.  find next 0.
            	  {
            	  if (fabs(cos(LMCConst::Pi()/2.*round((phi_shortTM01+LMCConst::Pi())/(LMCConst::Pi()/2.)))) < 1.e-4)  // find next zero.
            	    {
            	     dphi = LMCConst::Pi()/2.*round((phi_shortTM01+LMCConst::Pi())/(LMCConst::Pi()/2.)) - phi_shortTM01;
            	     phi_shortTM01 += dphi;
            	    }
            	  else  // try again.
            	    {
            	     dphi = LMCConst::Pi()/2.*round((phi_shortTM01+LMCConst::Pi()/2.)/(LMCConst::Pi()/2.)) - phi_shortTM01;
            	     phi_shortTM01 += dphi;
            	    }
            	  // continue toward polarizer, stop at electron and calculate field there:
            	  FieldFromShort += cos(phi_shortTM01 + 2.*LMCConst::Pi()*(CENTER_TO_SHORT + tPositionZ)/lambda_short);
            	  } // end short.
            }

            for (int i=0; i<nbounces; i++)  // polarizer-going wave, initially.  It already hit the polarizer but no one cares because the phase will just reset at the short.
            {
            	if (i%2==0)  // toward short.
            	  {
            	  if (fabs(cos(LMCConst::Pi()/2.*round((phi_polarizerTM01+LMCConst::Pi())/(LMCConst::Pi()/2.)))) < 1.e-4)  // find next zero.
            	    {
            	      dphi = LMCConst::Pi()/2.*round((phi_polarizerTM01+LMCConst::Pi())/(LMCConst::Pi()/2.)) - phi_polarizerTM01;
             	      phi_polarizerTM01 += dphi;
             	    }
            	  else
            	    {
            	      dphi = LMCConst::Pi()/2.*round((phi_polarizerTM01+LMCConst::Pi()/2.)/(LMCConst::Pi()/2.)) - phi_polarizerTM01;
             	      phi_polarizerTM01 += dphi;
            	    }
            	  // continue toward polarizer, stop at electron and calculate field there:
            	  FieldFromPolarizer += cos(phi_polarizerTM01 + 2.*LMCConst::Pi()*(CENTER_TO_SHORT + tPositionZ)/lambda_polarizer);
            	  } // end toward short loop.
            	else  // toward polarizer.
            	  {
            	  phi_polarizerTM01 += 2.*LMCConst::Pi()*(CENTER_TO_SHORT + CENTER_TO_ANTENNA)/lambda_polarizer + LMCConst::Pi();  // phase shift.
            	  // BACK UP toward short, stop at electron and calculate field there:
            	  FieldFromPolarizer += cos(phi_polarizerTM01 - 2.*LMCConst::Pi()*(CENTER_TO_ANTENNA - tPositionZ)/lambda_polarizer);
            	  }
            }

	    TM01FieldAfterBounces = (FieldFromShort + FieldFromPolarizer);

    	return TM01FieldAfterBounces;
    }


    double CyclotronRadiationExtractor::GetTM01FieldWithTerminator(KSParticle& anInitialParticle, KSParticle& aFinalParticle)
    {
      double dt = aFinalParticle.GetTime() - anInitialParticle.GetTime();
      double tCyclotronFrequency = aFinalParticle.GetCyclotronFrequency();
      double GroupVelocity = GetGroupVelocityTM01(aFinalParticle);
      double tVelocityZ = aFinalParticle.GetVelocity().GetZ();
      double tPositionZ = aFinalParticle.GetPosition().GetZ();
      double GammaZ = 1.0/sqrt(1.0-pow(tVelocityZ / GetGroupVelocityTM01(aFinalParticle),2.) );
      //      printf("tcyc is %g and GammaZ is %g and tVelocityZ is %g\n", tCyclotronFrequency, GammaZ, tVelocityZ);             
      double fprime_polarizer = tCyclotronFrequency*GammaZ*(1.-tVelocityZ/GroupVelocity);
      //      printf("GroupVelocity is %.10g and fprime_short is %.10g\n", GroupVelocity, fprime_short);                         
      double lambda_polarizer = GroupVelocity/fprime_polarizer;
      double FieldFromPolarizer=0.; // other doppler shift                                                               
      double TM01FieldAfterBounces = 0.;
      double phi_polarizerTM01 = 0.;

      phi_polarizerTM01 = 2.*LMCConst::Pi()*(2.*(CENTER_TO_ANTENNA-tPositionZ))
	/(GroupVelocity/fprime_polarizer);
      TM01FieldAfterBounces = cos(0.) + cos(phi_polarizerTM01);
      //printf("TE11FieldAfterOneBounce is %f\n", TE11FieldAfterOneBounce);                                              
      return TM01FieldAfterBounces;

 
    }


    double CyclotronRadiationExtractor::GetDampingFactorPhase1(KSParticle& anInitialParticle, KSParticle& aFinalParticle)
    {
        double TE01FieldFromShort = GetTE01FieldAfterOneBounce(anInitialParticle, aFinalParticle);
        double CouplingFactorTE01 = GetCouplingFactorTE01(aFinalParticle);
        double DampingFactorTE01 = CouplingFactorTE01*(1. - TE01FieldFromShort*TE01FieldFromShort);  // can be > 0 or < 0.

    	return DampingFactorTE01;
    }




    double CyclotronRadiationExtractor::GetDampingFactorPhase2(KSParticle& anInitialParticle, KSParticle& aFinalParticle)
    {
        double TE11FieldFromShort = GetTE11FieldAfterOneBounce(anInitialParticle, aFinalParticle);
        double TM01FieldWithTerminator = GetTM01FieldWithTerminator(anInitialParticle, aFinalParticle);
        double CouplingFactorTE11 = GetCouplingFactorTE11(aFinalParticle);
        double CouplingFactorTM01 = GetCouplingFactorTM01(aFinalParticle);

	//        double DampingFactorTE11 = CouplingFactorTE11*(1. - TE11FieldFromShort*TE11FieldFromShort);  // can be > 0 or < 0.
        double DampingFactorTM01 = CouplingFactorTM01*(1. - TM01FieldWithTerminator*TM01FieldWithTerminator);  // can be > 0 or < 0.
        double DampingFactor = DampingFactorTM01;

    	return DampingFactor;
    }


    void CyclotronRadiationExtractor::SetTrajectory( Kassiopeia::KSTrajectory* aTrajectory )
    {
        //this function is being run by KSRoot.cxx at initialization but presently
        //fTrajectory is not staying defined through the stepping.  Suspect binding problem.

    	fTrajectory = aTrajectory;

        return;
    }

    locust::Particle CyclotronRadiationExtractor::ExtractKassiopeiaParticle( KSParticle &anInitialParticle, KSParticle &aFinalParticle)
    {
        LMCThreeVector tPosition(aFinalParticle.GetPosition().Components());
        LMCThreeVector tVelocity(aFinalParticle.GetVelocity().Components());
        LMCThreeVector tMagneticField(aFinalParticle.GetMagneticField().Components());
        double tMass = aFinalParticle.GetMass();
        double tCharge = aFinalParticle.GetCharge();
        double tCyclotronFrequency = aFinalParticle.GetCyclotronFrequency();
        double tPitchAngle = aFinalParticle.GetPolarAngleToB();
    	double tTime = aFinalParticle.GetTime();


        locust::Particle aNewParticle;
        aNewParticle.SetPosition(tPosition.X(),tPosition.Y(),tPosition.Z());
        aNewParticle.SetVelocityVector(tVelocity.X(),tVelocity.Y(),tVelocity.Z());
        aNewParticle.SetMagneticFieldVector(tMagneticField.X(),tMagneticField.Y(),tMagneticField.Z());
        aNewParticle.SetMass(tMass);
        aNewParticle.SetCharge(tCharge);
        aNewParticle.SetTime(tTime);
        aNewParticle.SetCyclotronFrequency(2.*LMCConst::Pi()*tCyclotronFrequency);
        aNewParticle.SetKinematicProperties();

        if (fPitchAngle == -99.)  // first crossing of center
        {
        if (anInitialParticle.GetPosition().GetZ()/aFinalParticle.GetPosition().GetZ() < 0.)  // trap center
          {
          fPitchAngle = aFinalParticle.GetPolarAngleToB();
//      	printf("pitch angle is %f\n", fPitchAngle); getchar();

          }
        }
        aNewParticle.SetPitchAngle(fPitchAngle);


        return aNewParticle;

    }


    bool CyclotronRadiationExtractor::ExecutePostStepModification( KSParticle& anInitialParticle, KSParticle& aFinalParticle, KSParticleQueue& aQueue )
    {
        double DeltaE=0.;
        if(fP8Phase==1)
        {
            // adjust power with reflections.
            DeltaE = GetDampingFactorPhase1(anInitialParticle, aFinalParticle)*(aFinalParticle.GetKineticEnergy() - anInitialParticle.GetKineticEnergy());
            aFinalParticle.SetKineticEnergy((aFinalParticle.GetKineticEnergy() - DeltaE));
        }
        if(fP8Phase==2)
        {
            // adjust power with reflections.
            DeltaE = GetDampingFactorPhase2(anInitialParticle, aFinalParticle)*(aFinalParticle.GetKineticEnergy() - anInitialParticle.GetKineticEnergy());
            //printf("poststep says DeltaE is %g\n", DeltaE);
            aFinalParticle.SetKineticEnergy((aFinalParticle.GetKineticEnergy() - DeltaE));
        }

        if (!fDoneWithSignalGeneration)  // if Locust is still acquiring voltages.
        {

        if (t_old == 0.) fPitchAngle = -99.;  // new electron needs central pitch angle reset.
    	double t_poststep = aFinalParticle.GetTime();
        fNewParticleHistory.push_back(ExtractKassiopeiaParticle(anInitialParticle, aFinalParticle));

        if (t_poststep - t_old >= fDigitizerTimeStep) //take a digitizer sample every 5e-10s
        {
            std::unique_lock< std::mutex >tLock( fMutexDigitizer, std::defer_lock );  // lock access to mutex before writing to globals.
            tLock.lock();

            int tHistoryMaxSize;

            //Dont want to check .back() of history if it is empty! -> Segfault
            if(fParticleHistory.size() && (fNewParticleHistory.back().GetTime() < fParticleHistory.back().GetTime()))
            {
//                printf("New Particle!, t_old is %g\n", t_old); getchar();
                t_poststep = 0.;
                fParticleHistory.clear();
            }

            //Phase I or II Setup: Put only last particle in fParticleHistory. Use interpolated value for the particle
            if((fP8Phase==2) || (fP8Phase==1))
            {
                // interpolate particle state.  Have to pull trajectory out of toolbox due to binding problem in SetTrajectory above.
                KSParticle tParticleCopy = aFinalParticle;
                katrin::KToolbox::GetInstance().Get< Kassiopeia::KSTrajectory  >( "root_trajectory" )->GetInterpolatedParticleState(t_old + fDigitizerTimeStep, tParticleCopy);
                fParticleHistory.push_back(ExtractKassiopeiaParticle(anInitialParticle, tParticleCopy));

                tHistoryMaxSize = 5;

            }
            else
            {
                //Put in new entries in global ParticleHistory
                fParticleHistory.insert(fParticleHistory.end(),fNewParticleHistory.begin(),fNewParticleHistory.end());

                for(int i=fParticleHistory.size()-fNewParticleHistory.size()-1;i<fParticleHistory.size()-1;i++)
                {
                    fParticleHistory[i].SetSpline(fParticleHistory[i+1]);
                }

                tHistoryMaxSize = 5000;

            }

            fNewParticleHistory.clear();

            //Purge fParticleHistory of overly old entries	    
	    while(t_poststep-fParticleHistory.front().GetTime()>1e-7 || fParticleHistory.size() > tHistoryMaxSize)
	      {
	      fParticleHistory.pop_front();
	      }
	     //	    printf("done purging\n");
	    
            tLock.unlock();
            fDigitizerCondition.notify_one();  // notify Locust after writing.

        }
        } // fDoneWithSignalGeneration

        return true;
    }


    void CyclotronRadiationExtractor::InitializeComponent()
    {
    }

    void CyclotronRadiationExtractor::DeinitializeComponent()
    {
    }

    void CyclotronRadiationExtractor::PullDeupdateComponent()
    {
    }
    void CyclotronRadiationExtractor::PushDeupdateComponent()
    {
    }

} /* namespace locust */
