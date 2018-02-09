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
         double SpeedOfLight = 2.99792458e8; // m/s
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
    	return coupling;
    }

    double CyclotronRadiationExtractor::GetCouplingFactorTM01(KSParticle& aFinalParticle)
    {
    	double kc = 2.405/0.00502920;
    	double x = aFinalParticle.GetPosition().GetX();
    	double y = aFinalParticle.GetPosition().GetY();
    	double r = sqrt(x*x + y*y);
    	double coupling =   146876.5/168.2 * 2./LMCConst::Pi() * 4./(2.*LMCConst::Pi()) / kc * j1(kc*r);
    	return coupling;
    }

    double CyclotronRadiationExtractor::GetCouplingFactorTE01(KSParticle& aFinalParticle)  // Phase 1
    {
    	double dim1_wr42 = 10.668e-3; // a in m
    	double x = aFinalParticle.GetPosition().GetX() + dim1_wr42/2.;
    	double coupling = 0.63*sin(LMCConst::Pi()*x/dim1_wr42);  // avg over cyclotron orbit.
    	return coupling;
    }


    double CyclotronRadiationExtractor::GetTE01FieldAfterOneBounce(KSParticle& anInitialParticle, KSParticle& aFinalParticle)
    {
        double fcyc = aFinalParticle.GetCyclotronFrequency();
        double GroupVelocity = GetGroupVelocityTE01(aFinalParticle);
    	double zvelocity = aFinalParticle.GetVelocity().GetZ();
    	double zPosition = aFinalParticle.GetPosition().GetZ();
        double GammaZ = 1.0/pow(1.0-pow(zvelocity/GetGroupVelocityTE01(aFinalParticle),2.),0.5);

    	double fprime_short = fcyc*GammaZ*(1.+zvelocity/GroupVelocity);
    	double phi_short = 2.*LMCConst::Pi()*2.*(zPosition+CENTER_TO_SHORT)/(GroupVelocity/fprime_short);
//        double FieldFromShort = cos(phi_short);  // no resonant enhancement.
        double FieldFromShort = cos(0.) + cos(phi_short); // yes resonant enhancement.

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

        phi_shortTE11 = 2.*LMCConst::Pi()*2.*(zPosition+CENTER_TO_SHORT)/(GroupVelocity/fprime_short);
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
    	int nbounces = 20;
        double time_decay = 1.;
        double reflection_coefficient = 1.0;
        double phi_shortTM01 = 0.;
        double phi_polarizerTM01 = 0.;

        //printf("TM01 l1 is %g and l2 is %g\n", GroupVelocity/fprime_short, GroupVelocity/fprime_polarizer); getchar();

//        if ((phi_shortTM01[0] == 0.)||(0==0))  // if the event has just started, or always.
//        {
            phi_shortTM01 = 2.* LMCConst::Pi() *2.*(tPositionZ+CENTER_TO_SHORT)/lambda_short;  // starting phi after 0th bounce.
            phi_polarizerTM01 = 2.*LMCConst::Pi()*2.*(CENTER_TO_ANTENNA - tPositionZ)/lambda_polarizer + LMCConst::Pi();  // starting phi after 0th bounce.

//   	       	printf("phi_shortTM01[0] is %.10g and LMCConst::Pi() is %.10g and z is %.10g and lambda is %.10g\n", phi_shortTM01[0], LMCConst::Pi(), tPositionZ, lambda_short);

            FieldFromShort = cos(0.) + 1./1.4*reflection_coefficient*cos(phi_shortTM01); // starting field, after 0th bounce.
            FieldFromPolarizer = 1./1.4*reflection_coefficient*cos(phi_polarizerTM01); // starting field, after 0th bounce.


            for (int i=0; i<nbounces; i++)  // short-going wave, initially.
            {
	      //	      time_decay = exp(-(double)i*2./(double)nbounces);
                if (i%2==0)
                {
                    phi_shortTM01 += 2.*LMCConst::Pi()*2.*(CENTER_TO_ANTENNA - tPositionZ)/lambda_short + LMCConst::Pi();  // phase shift PI
                }
                else
                {
                    phi_shortTM01 += 2.*LMCConst::Pi()*2.*(tPositionZ + CENTER_TO_SHORT)/lambda_short;
                }

                FieldFromShort += time_decay*cos(phi_shortTM01); // field adds after each bounce.
            }


            for (int i=0; i<nbounces; i++)  // polarizer-going wave, initially.
            {
	      //	      time_decay = exp(-(double)i*2./(double)nbounces);
                if (i%2==0)
                {
                    phi_polarizerTM01 += 2.*LMCConst::Pi()*2.*(CENTER_TO_SHORT + tPositionZ)/lambda_polarizer;
                }
                else
                {
                    phi_polarizerTM01 += phi_polarizerTM01 + 2.*LMCConst::Pi()*2.*(CENTER_TO_ANTENNA - tPositionZ)/lambda_polarizer + LMCConst::Pi();  // phase shift PI.
                }
                FieldFromPolarizer += time_decay*cos(phi_polarizerTM01);
            }
//        } // phi_shortTM01 == 0.  This loop defines the initial standing wave on the first tracking step.

        TM01FieldAfterBounces = FieldFromShort + FieldFromPolarizer;

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
        double TM01FieldAfterBounces = GetTM01FieldAfterBounces(anInitialParticle, aFinalParticle);
        double CouplingFactorTE11 = GetCouplingFactorTE11(aFinalParticle);
        double CouplingFactorTM01 = GetCouplingFactorTM01(aFinalParticle);

        double DampingFactorTE11 = CouplingFactorTE11*(1. - TE11FieldFromShort*TE11FieldFromShort);  // can be > 0 or < 0.
        double DampingFactorTM01 = CouplingFactorTM01*(1. - TM01FieldAfterBounces*TM01FieldAfterBounces);  // can be > 0 or < 0.
        double DampingFactor = DampingFactorTM01 + DampingFactorTE11;

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
        KGeoBag::KThreeVector tPosition = aFinalParticle.GetPosition();
        KGeoBag::KThreeVector tVelocity = aFinalParticle.GetVelocity();
        KGeoBag::KThreeVector tMagneticField = aFinalParticle.GetMagneticField();
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
          }
        }
        aNewParticle.SetPitchAngle(fPitchAngle);


        return aNewParticle;

    }


    bool CyclotronRadiationExtractor::ExecutePostStepModification( KSParticle& anInitialParticle, KSParticle& aFinalParticle, KSParticleQueue& aQueue )
    {
//      printf("fcyc before coupling is %.9g and Bz is %.10g\n\n", aFinalParticle.GetCyclotronFrequency(), aFinalParticle.GetMagneticField().GetZ());

        //printf("pre step kinetic energy - 4.84338e-15 is %g\n", anInitialParticle.GetKineticEnergy()- 4.84338e-15); //getchar();
        //printf("post step kinetic energy - 4.84338e-15 is %g\n", aFinalParticle.GetKineticEnergy()- 4.84338e-15); //getchar();

        double DeltaE=0.;

        if(fP8Phase==1)
        {
            // adjust power with reflections.
            DeltaE = GetDampingFactorPhase1(anInitialParticle, aFinalParticle)*(aFinalParticle.GetKineticEnergy() - anInitialParticle.GetKineticEnergy());
            //printf("poststep says DeltaE is %g\n", DeltaE);
            aFinalParticle.SetKineticEnergy((aFinalParticle.GetKineticEnergy() - DeltaE));
        }
        if(fP8Phase==2)
        {
            // adjust power with reflections.
            DeltaE = GetDampingFactorPhase2(anInitialParticle, aFinalParticle)*(aFinalParticle.GetKineticEnergy() - anInitialParticle.GetKineticEnergy());
            //printf("poststep says DeltaE is %g\n", DeltaE);
            aFinalParticle.SetKineticEnergy((aFinalParticle.GetKineticEnergy() - DeltaE));
        }

    	double t_poststep = aFinalParticle.GetTime();
        fNewParticleHistory.push_back(ExtractKassiopeiaParticle(anInitialParticle, aFinalParticle));

        if (t_poststep - t_old >= fDigitizerTimeStep) //take a digitizer sample every 5e-10s
        {
            std::unique_lock< std::mutex >tLock( fMutexDigitizer, std::defer_lock );  // lock access to mutex before writing to globals.
            tLock.lock();

            int tHistoryMaxSize;

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
                fParticleHistory.pop_front();

            tLock.unlock();
            fDigitizerCondition.notify_one();  // notify Locust after writing.

        }

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
