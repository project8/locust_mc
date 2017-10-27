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

    CyclotronRadiationExtractor::CyclotronRadiationExtractor()
    {
    }

    CyclotronRadiationExtractor::CyclotronRadiationExtractor( const CyclotronRadiationExtractor& aOrig ) 
    {
    }

    CyclotronRadiationExtractor::~CyclotronRadiationExtractor()
    {
    }

    CyclotronRadiationExtractor* CyclotronRadiationExtractor::Clone() const
    {
        return new CyclotronRadiationExtractor( *this );
    }


    bool CyclotronRadiationExtractor::ExecutePreStepModification( KSParticle& anInitialParticle, KSParticleQueue& aQueue )
    {
    	return true;
    }

    double CyclotronRadiationExtractor::GetGroupVelocityTM01(KSParticle& aFinalParticle)
    {
        double SpeedOfLight = KConst::C();
        double CutOffFrequency = 2. * KConst::Pi() * SpeedOfLight * 2.405 / 2. / KConst::Pi() / 0.00502920; // rad/s
        double fcyc = aFinalParticle.GetCyclotronFrequency();
        double GroupVelocity = SpeedOfLight * sqrt( 1. - pow(CutOffFrequency/(2.*KConst::Pi()*fcyc), 2.) );
//        printf("GroupVelocity is %g\n", GroupVelocity); getchar();
    	return GroupVelocity;
    }
    double CyclotronRadiationExtractor::GetGroupVelocityTE11(KSParticle& aFinalParticle)
    {
        double SpeedOfLight = KConst::C(); // m/s
        double CutOffFrequency = 2. * KConst::Pi() * SpeedOfLight * 1.841 / 2. / KConst::Pi() / 0.00502920; // rad/s
        double fcyc = aFinalParticle.GetCyclotronFrequency();
        double GroupVelocity = SpeedOfLight * sqrt( 1. - pow(CutOffFrequency/(2.*KConst::Pi()*fcyc), 2.) );
        //("GroupVelocity is %g\n", GroupVelocity); getchar();
    	return GroupVelocity;
    }


    double CyclotronRadiationExtractor::GetCouplingFactorTE11(KSParticle& aFinalParticle)
    {
    	double kc = 1.841/0.00502920;
    	double x = aFinalParticle.GetPosition().GetX();
    	double y = aFinalParticle.GetPosition().GetY();
    	double dx = 0.;
    	double dy = 0.;
        //double kassgcpx = aFinalParticle.GetGuidingCenterPosition().GetX();
        //double kassgcpy = aFinalParticle.GetGuidingCenterPosition().GetY();
        //double vx = aFinalParticle.GetVelocity().GetX();
        //double vy = aFinalParticle.GetVelocity().GetY();
        //double v = aFinalParticle.GetTransVelocity();
        //double cyc_radius = aFinalParticle.GetTransVelocity()/(2.*PI*aFinalParticle.GetCyclotronFrequency());
        //if (vx/vy > 0.) dx = cyc_radius*vx/v;
        //else dx = -cyc_radius*vx/v;
        //if (vx/vy < 0.) dy = cyc_radius*vy/v;
        //else dy = -cyc_radius*vy/v;

        //x = x + dx;
        //y = y + dy;
        //printf("cyc_radius is %f\n", cyc_radius);
        //printf("x is %f and y is %f\n", x, y);

        //printf("gcpx = %f and gcpy = %f\n", gcpx, gcpy);
        //printf("kass gcpx = %f and kass gcpy = %f\n", kassgcpx, kassgcpy); getchar();
        //printf("dx is %f and dy is %f\n", dx, dy);
    	double r = sqrt( x*x + y*y);
    	double coupling = 119116./168.2 * 2./KConst::Pi() * 4./(2.*KConst::Pi()) / kc/2. * ( (j0(kc*r) - jn(2,kc*r)) +
    			(j0(kc*r) + jn(2, kc*r)) );
        //printf("TE11 coupling at r=%f is %f\n", r, coupling); //getchar();
    	return coupling;
    }

    double CyclotronRadiationExtractor::GetCouplingFactorTM01(KSParticle& aFinalParticle)
    {
    	double kc = 2.405/0.00502920;
    	double x = aFinalParticle.GetPosition().GetX();
    	double y = aFinalParticle.GetPosition().GetY();
    	double r = sqrt(x*x + y*y);
    	double coupling =   146876.5/168.2 * 2./KConst::Pi() * 4./(2.*KConst::Pi()) / kc * j1(kc*r);

    	printf("j1(kc*r) is %f\n", gsl_sf_bessel_zero_Jnu(0,1)); getchar();

    	//printf("tm01 coupling at r=%f is %f\n", r, coupling); //getchar();
        //printf("guiding center z is %f\n", aFinalParticle.GetGuidingCenterPosition().GetZ());
    	return coupling;
    }



    double CyclotronRadiationExtractor::GetTE11FieldAfterOneBounce(KSParticle& anInitialParticle, KSParticle& aFinalParticle)
    {
    	double dt = aFinalParticle.GetTime() - anInitialParticle.GetTime();
        double fcyc = aFinalParticle.GetCyclotronFrequency();
        double GroupVelocity = GetGroupVelocityTE11(aFinalParticle);

    	double zvelocity = aFinalParticle.GetVelocity().GetZ();
    	double zposition = aFinalParticle.GetPosition().GetZ();
        double GammaZ = 1.0/pow(1.0-pow(zvelocity/GroupVelocity,2.),0.5);
    	double fprime_short = fcyc*GammaZ*(1.+zvelocity/GroupVelocity);
    	double TE11FieldAfterOneBounce = 0.;

    	if ((phi_shortTE11==0.)||(0==0))  // after each step.
        {
            phi_shortTE11 = 2.*KConst::Pi()*2.*(zposition+CENTER_TO_SHORT)/(GroupVelocity/fprime_short);
            TE11FieldAfterOneBounce = cos(0.) + cos(phi_shortTE11);
            //printf("TE11FieldAfterOneBounce is %f\n", TE11FieldAfterOneBounce);
        }
    	//else
        //{
        //    phi_shortTE11 += 2.*PI*fprime_short*dt;
        //    TE11FieldAfterOneBounce = cos(0.) + cos(phi_shortTE11);
    	//  
        //}

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
        double tPI = KConst::Pi();
    	int nbounces = 10;

        //printf("TM01 l1 is %g and l2 is %g\n", GroupVelocity/fprime_short, GroupVelocity/fprime_polarizer); getchar();

        if ((phi_shortTM01[0] == 0.)||(0==0))  // if the event has just started, or always.
        {
            phi_shortTM01[0] = 2.* tPI *2.*(tPositionZ+CENTER_TO_SHORT)/lambda_short;  // starting phi after 0th bounce.
            phi_polarizerTM01[0] = 2.*tPI*2.*(CENTER_TO_ANTENNA - tPositionZ)/lambda_polarizer + tPI;  // starting phi after 0th bounce.

//   	       	printf("phi_shortTM01[0] is %.10g and tPI is %.10g and z is %.10g and lambda is %.10g\n", phi_shortTM01[0], tPI, tPositionZ, lambda_short);

            FieldFromShort = cos(0.) + cos(phi_shortTM01[0]); // starting field, after 0th bounce.
            FieldFromPolarizer = cos(0.) + cos(phi_polarizerTM01[0]); // starting field, after 0th bounce.
//   	       	printf("starting TM01FieldFromShort is %g\n", FieldFromShort);


            for (int i=0; i<nbounces; i++)  // short-going wave, initially.
            {
                if (i%2==0)
                {
                    phi_shortTM01[i+1] = phi_shortTM01[i] + 2.*tPI*2.*(CENTER_TO_ANTENNA - tPositionZ)/lambda_short + tPI;  // phase shift
                }
                else
                {
                    phi_shortTM01[i+1] = phi_shortTM01[i] + 2.*tPI*2.*(tPositionZ + CENTER_TO_SHORT)/lambda_short;
                }
                FieldFromShort += cos(phi_shortTM01[i+1]); // field adds after each bounce.
//       	       	printf("TM01FieldFromShort is %g\n", FieldFromShort);

            }


            for (int i=0; i<nbounces; i++)  // polarizer-going wave, initially.
            {
                if (i%2==0)
                {
                    phi_polarizerTM01[i+1] = phi_polarizerTM01[i] + 2.*tPI*2.*(CENTER_TO_SHORT + tPositionZ)/lambda_polarizer;
                }
                else
                {
                    phi_polarizerTM01[i+1] = phi_polarizerTM01[i] + 2.*tPI*2.*(CENTER_TO_ANTENNA - tPositionZ)/lambda_polarizer + tPI;  // phase shift.
                }
                FieldFromPolarizer += cos(phi_polarizerTM01[i+1]);
            }
        } // phi_shortTM01 == 0.  This loop defines the initial standing wave on the first tracking step.


//        printf("fieldfromshort is %.10g\n", FieldFromShort);
//        printf("fieldfrompolarizer is %.10g\n", FieldFromPolarizer);
//        getchar();

        TM01FieldAfterBounces = FieldFromShort + FieldFromPolarizer;

        //printf("x y z is %f %f %f\n", aFinalParticle.GetPosition().GetX(), aFinalParticle.GetPosition().GetY(), tPositionZ);
        //printf("lambda_short is %g and lambda_polarizer is %g\n", lambda_short, lambda_polarizer);
        //printf("TM01FieldAfterBounces is %g\n", TM01FieldAfterBounces); //getchar();

        //printf("TM01FieldAfterBounces is %f\n", TM01FieldAfterBounces);

    	return TM01FieldAfterBounces;
    }



    double CyclotronRadiationExtractor::GetDampingFactor(KSParticle& anInitialParticle, KSParticle& aFinalParticle)
    {
        //double fcyc = aFinalParticle.GetCyclotronFrequency();
        //double GroupVelocity = GetGroupVelocityTE11(aFinalParticle);
        //double zvelocity = aFinalParticle.GetVelocity().GetZ();
        //double zposition = aFinalParticle.GetPosition().GetZ();
        //double GammaZ = 1.0/pow(1.0-pow(zvelocity/GetGroupVelocityTE11(aFinalParticle),2.),0.5);

        //double fprime_short = fcyc*GammaZ*(1.+zvelocity/GroupVelocity);

        //double phi_short = 2.*PI*2.*(zposition+CENTER_TO_SHORT)/(GroupVelocity/fprime_short);
        //double TE11FieldFromShort = cos(0.) + cos(phi_short); // resonant enhancement.
        double TE11FieldFromShort = GetTE11FieldAfterOneBounce(anInitialParticle, aFinalParticle);
        double TM01FieldAfterBounces = GetTM01FieldAfterBounces(anInitialParticle, aFinalParticle);
        double CouplingFactorTE11 = GetCouplingFactorTE11(aFinalParticle);
        double CouplingFactorTM01 = GetCouplingFactorTM01(aFinalParticle);

//        printf("TE11FieldFromShort is %.10g and TM01FieldAfterBounces is %.10g and CouplingTE11 is %.10g and CouplingTM01 is %.10g\n",
//        		TE11FieldFromShort, TM01FieldAfterBounces, CouplingFactorTE11, CouplingFactorTM01); getchar();

        double DampingFactorTE11 = CouplingFactorTE11*(1. - TE11FieldFromShort*TE11FieldFromShort);  // can be > 0 or < 0.
        double DampingFactorTM01 = CouplingFactorTM01*(1. - TM01FieldAfterBounces*TM01FieldAfterBounces);  // can be > 0 or < 0.
        double DampingFactor = DampingFactorTM01 + DampingFactorTE11;
		//double DampingFactor = DampingFactorTE11;

        //if (fabs(DampingFactor)>0.)
        //{
        //    printf("elapsed time is %g\n", aFinalParticle.GetTime());
        //    printf("x, y, z position is %f %f %f\n", aFinalParticle.GetPosition().GetX(), aFinalParticle.GetPosition().GetY(), aFinalParticle.GetPosition().GetZ());
        //    printf("couplingTM01 is %f\n", CouplingFactorTM01);
        //    printf("damping factor total = %f, DampingFactorTE11 is %g and DampingFactorTM01 is %g\n", DampingFactor, DampingFactorTE11, DampingFactorTM01);
        //    getchar();
        //}
                    

    	return DampingFactor;
    }


    void CyclotronRadiationExtractor::SetTrajectory( Kassiopeia::KSTrajectory* aTrajectory )
    {
        //this function is being run by KSRoot.cxx at initialization but presently
        //fTrajectory is not staying defined through the stepping.  Suspect binding problem.

    	fTrajectory = aTrajectory;

        return;
    }

    locust::Particle CyclotronRadiationExtractor::ExtractKassiopeiaParticle( KSParticle &aFinalParticle)
    {
        KGeoBag::KThreeVector tPosition = aFinalParticle.GetPosition();
        KGeoBag::KThreeVector tVelocity = aFinalParticle.GetVelocity();
        KGeoBag::KThreeVector tMagneticField = aFinalParticle.GetMagneticField();
        double tMass = aFinalParticle.GetMass();
        double tCharge = aFinalParticle.GetCharge();
        double tCyclotronFrequency = aFinalParticle.GetCyclotronFrequency();
    	double tTime = aFinalParticle.GetTime();

        //printf("%e \n",fcyc);

        locust::Particle aNewParticle;
        aNewParticle.SetPosition(tPosition.X(),tPosition.Y(),tPosition.Z());
        aNewParticle.SetVelocityVector(tVelocity.X(),tVelocity.Y(),tVelocity.Z());
        aNewParticle.SetMagneticFieldVector(tMagneticField.X(),tMagneticField.Y(),tMagneticField.Z());
        aNewParticle.SetMass(tMass);
        aNewParticle.SetCharge(tCharge);
        aNewParticle.SetTime(tTime);
        aNewParticle.SetCyclotronFrequency(2.*KConst::Pi()*tCyclotronFrequency);
        aNewParticle.SetKinematicProperties();
        
        return aNewParticle;

    }


    bool CyclotronRadiationExtractor::ExecutePostStepModification( KSParticle& anInitialParticle, KSParticle& aFinalParticle, KSParticleQueue& aQueue )
    {
//      printf("fcyc before coupling is %.9g and vz is %.10g\n\n", aFinalParticle.GetCyclotronFre\
quency(), aFinalParticle.GetVelocity().GetZ()); getchar();

        //printf("pre step kinetic energy - 4.84338e-15 is %g\n", anInitialParticle.GetKineticEnergy()- 4.84338e-15); //getchar();
        //printf("post step kinetic energy - 4.84338e-15 is %g\n", aFinalParticle.GetKineticEnergy()- 4.84338e-15); //getchar();

        double DeltaE=0.;
        if(fPhaseIISimulation)
        {
            // adjust power with reflections.
            DeltaE = GetDampingFactor(anInitialParticle, aFinalParticle)*(aFinalParticle.GetKineticEnergy() - anInitialParticle.GetKineticEnergy());
            //printf("poststep says DeltaE is %g\n", DeltaE);
            aFinalParticle.SetKineticEnergy((aFinalParticle.GetKineticEnergy() - DeltaE));
        }


//        printf("z is %f and DeltaE is %g and post fix kinetic energy is %g and fcyc is %.9g\n", aFinalParticle.GetPosition().Z(), DeltaE, aFinalParticle.GetKineticEnergy() - 4.84338e-15, aFinalParticle.GetCyclotronFrequency()); getchar();

    	t_poststep = aFinalParticle.GetTime();

        fNewParticleHistory.push_back(ExtractKassiopeiaParticle(aFinalParticle));

        if (t_poststep - t_old >= fDigitizerTimeStep) //take a digitizer sample every 5e-10s
        {
            std::unique_lock< std::mutex >tLock( fMutexDigitizer, std::defer_lock );  // lock access to mutex before writing to globals.
            tLock.lock();

//            t_poststep = t_old + fDigitizerTimeStep;


            int tHistoryMaxSize;

            //Phase II Setup: Put only last particle in fParticleHistory. Use interpolated value for the particle
            if(fPhaseIISimulation)
            {
                // interpolate particle state.  Have to pull trajectory out of toolbox due to binding problem in SetTrajectory above.
                KSParticle tParticleCopy = aFinalParticle;
                katrin::KToolbox::GetInstance().Get< Kassiopeia::KSTrajectory  >( "root_trajectory" )->GetInterpolatedParticleState(t_old+fDigitizerTimeStep, tParticleCopy);
                fParticleHistory.push_back(ExtractKassiopeiaParticle(tParticleCopy));
//                printf("t_old is %g and t_poststep is %g\n", t_old, t_poststep);
//                printf("\n\n\nAbout to digitize:  z is %g and DeltaE is %g and post fix kinetic energy is %g and fcyc is %.9g\n",
//                		tParticleCopy.GetPosition().Z(), DeltaE, tParticleCopy.GetKineticEnergy() - 4.84338e-15, tParticleCopy.GetCyclotronFrequency()); getchar();

                tHistoryMaxSize = 5;

            }
            else
            {
                //Put in new entries in global ParticleHistory
                fParticleHistory.insert(fParticleHistory.end(),fNewParticleHistory.begin(),fNewParticleHistory.end());
                //Set Spline coefficients -
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

//            t_old = t_poststep;  // get ready to look for next sample.


             //printf("de is %g and dt is %g and LarmorPower is %g\n", de, dt, LarmorPower);
          	 //printf("Kassiopeia says:  tick has happened; continuous time is %g and zvelocity is %f\n", t_poststep, zvelocity);
          	 //printf("Kassiopeia says:  fcyc is %g\n", fcyc);
          	 //printf("  initial particle momentum: %g\n", anInitialParticle.GetMomentum().Z());
             //printf("Mass is %g\n", anInitialParticle.GetMass());
             //printf("k.e. = %g eV\n", anInitialParticle.GetKineticEnergy_eV());
             //printf("Lorentz factor is %f\n", anInitialParticle.GetLorentzFactor());
             //printf("1/(sqrt(1-v^2/c^2) is %f\n", 1.0/pow(1.0-pow(anInitialParticle.GetSpeed()/2.99792e8,2.),0.5));
             //printf("FinalParticle().IsActive() is %d\n", aFinalParticle.IsActive());
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
