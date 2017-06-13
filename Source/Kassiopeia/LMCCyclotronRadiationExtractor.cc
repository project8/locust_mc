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
        double SpeedOfLight = 2.99792458e8; // m/s
        double CutOffFrequency = 2. * PI * SpeedOfLight * 2.405 / 2. / PI / 0.00502920; // rad/s
        double fcyc = aFinalParticle.GetCyclotronFrequency();
        double GroupVelocity = SpeedOfLight * pow( 1. - pow(CutOffFrequency/(2.*PI*fcyc), 2.) , 0.5);
//        printf("GroupVelocity is %g\n", GroupVelocity); getchar();
    	return GroupVelocity;
    }

    double CyclotronRadiationExtractor::GetGroupVelocityTE11(KSParticle& aFinalParticle)
    {
        double SpeedOfLight = 2.99792458e8; // m/s
        double CutOffFrequency = 2. * PI * SpeedOfLight * 1.841 / 2. / PI / 0.00502920; // rad/s
        double fcyc = aFinalParticle.GetCyclotronFrequency();
        double GroupVelocity = SpeedOfLight * pow( 1. - pow(CutOffFrequency/(2.*PI*fcyc), 2.) , 0.5);
//        printf("GroupVelocity is %g\n", GroupVelocity); getchar();
    	return GroupVelocity;
    }


    double CyclotronRadiationExtractor::GetCouplingFactorTE11(KSParticle& aFinalParticle)
    {
    	double kc = 1.841/0.00502920;
    	double x = aFinalParticle.GetPosition().GetX();
    	double y = aFinalParticle.GetPosition().GetY();
    	double dx = 0.;
    	double dy = 0.;
//    	double kassgcpx = aFinalParticle.GetGuidingCenterPosition().GetX();
//    	double kassgcpy = aFinalParticle.GetGuidingCenterPosition().GetY();
//    	double vx = aFinalParticle.GetVelocity().GetX();
//    	double vy = aFinalParticle.GetVelocity().GetY();
 //   	double v = aFinalParticle.GetTransVelocity();
 //   	double cyc_radius = aFinalParticle.GetTransVelocity()/(2.*PI*aFinalParticle.GetCyclotronFrequency());
 //   	if (vx/vy > 0.) dx = cyc_radius*vx/v;
 //     	      else dx = -cyc_radius*vx/v;
 //      	if (vx/vy < 0.) dy = cyc_radius*vy/v;
 //         	  else dy = -cyc_radius*vy/v;

//    	x = x + dx;
//    	y = y + dy;
//    	printf("cyc_radius is %f\n", cyc_radius);
//    	printf("x is %f and y is %f\n", x, y);

//    	printf("gcpx = %f and gcpy = %f\n", gcpx, gcpy);
//    	printf("kass gcpx = %f and kass gcpy = %f\n", kassgcpx, kassgcpy); getchar();
//    	printf("dx is %f and dy is %f\n", dx, dy);
    	double r = pow(x*x+y*y,0.5);
    	double coupling = 119116./168.2 * 2./PI * 4./(2.*PI) / kc/2. * ( (j0(kc*r) - jn(2,kc*r)) +
    			(j0(kc*r) + jn(2, kc*r)) );
 //   	printf("TE11 coupling at r=%f is %f\n", r, coupling); //getchar();
    	return coupling;
    }

    double CyclotronRadiationExtractor::GetCouplingFactorTM01(KSParticle& aFinalParticle)
    {
    	double kc = 2.405/0.00502920;
    	double x = aFinalParticle.GetPosition().GetX();
    	double y = aFinalParticle.GetPosition().GetY();
    	double r = pow(x*x+y*y,0.5);
    	double coupling =   146876.5/168.2 * 2./PI * 4./(2.*PI) / kc * j1(kc*r);
//    	printf("tm01 coupling at r=%f is %f\n", r, coupling); //getchar();
//    	printf("guiding center z is %f\n", aFinalParticle.GetGuidingCenterPosition().GetZ());
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
      	  phi_shortTE11 = 2.*PI*2.*(zposition+CENTER_TO_SHORT)/(GroupVelocity/fprime_short);
      	  TE11FieldAfterOneBounce = cos(0.) + cos(phi_shortTE11);
//      	  printf("TE11FieldAfterOneBounce is %f\n", TE11FieldAfterOneBounce);
    	  }
    	else
    	  {
    		/*
    	  phi_shortTE11 += 2.*PI*fprime_short*dt;
    	  TE11FieldAfterOneBounce = cos(0.) + cos(phi_shortTE11);
    	  */
    	  }

    	return TE11FieldAfterOneBounce;
    }

    double CyclotronRadiationExtractor::GetTM01FieldAfterBounces(KSParticle& anInitialParticle, KSParticle& aFinalParticle)
    {

    	double dt = aFinalParticle.GetTime() - anInitialParticle.GetTime();
        double fcyc = aFinalParticle.GetCyclotronFrequency();
        double GroupVelocity = GetGroupVelocityTM01(aFinalParticle);
//       	double zvelocity = aFinalParticle.GetVelocity().GetZ();
        double zposition = aFinalParticle.GetPosition().GetZ();
        double GammaZ = 1.0/pow(1.0-pow(zvelocity/GetGroupVelocityTM01(aFinalParticle),2.),0.5);

    	double fprime_short = fcyc*GammaZ*(1.+zvelocity/GroupVelocity);
    	double fprime_polarizer = fcyc*GammaZ*(1.-zvelocity/GroupVelocity);
    	double lambda_short = GroupVelocity/fprime_short;
    	double lambda_polarizer = GroupVelocity/fprime_polarizer;
    	double FieldFromShort=0.;  // first doppler shift
    	double FieldFromPolarizer=0.; // other doppler shift
    	double TM01FieldAfterBounces = 0.;
    	int nbounces = 100;

//    	printf("TM01 l1 is %g and l2 is %g\n", GroupVelocity/fprime_short, GroupVelocity/fprime_polarizer); getchar();

        if ((phi_shortTM01[0] == 0.)||(0==0))  // if the event has just started, or always.
        {
    	phi_shortTM01[0] = 2.*PI*2.*(zposition+CENTER_TO_SHORT)/lambda_short;  // starting phi after 0th bounce.
    	phi_polarizerTM01[0] = 2.*PI*2.*(CENTER_TO_ANTENNA - zposition)/lambda_polarizer + PI;  // starting phi after 0th bounce.

    	FieldFromShort = cos(0.) + cos(phi_shortTM01[0]); // starting field, after 0th bounce.
    	FieldFromPolarizer = cos(0.) + cos(phi_polarizerTM01[0]); // starting field, after 0th bounce.

    	for (int i=0; i<nbounces; i++)  // short-going wave, initially.
    	{
    	if (i%2==0)
    	  {
    	  phi_shortTM01[i+1] = phi_shortTM01[i] + 2.*PI*2.*(CENTER_TO_ANTENNA - zposition)/lambda_short + PI;  // phase shift
    	  }
    	else
    	  {
    	  phi_shortTM01[i+1] = phi_shortTM01[i] + 2.*PI*2.*(zposition + CENTER_TO_SHORT)/lambda_short;
    	  }
  	    FieldFromShort += cos(phi_shortTM01[i+1]); // field adds after each bounce.
    	}


    	for (int i=0; i<nbounces; i++)  // polarizer-going wave, initially.
    	{
    	if (i%2==0)
    	  {
    	  phi_polarizerTM01[i+1] = phi_polarizerTM01[i] + 2.*PI*2.*(CENTER_TO_SHORT + zposition)/lambda_polarizer;
    	  }
    	else
    	  {
    	  phi_polarizerTM01[i+1] = phi_polarizerTM01[i] + 2.*PI*2.*(CENTER_TO_ANTENNA - zposition)/lambda_polarizer + PI;  // phase shift.
    	  }
  	    FieldFromPolarizer += cos(phi_polarizerTM01[i+1]);
    	}
        } // phi_shortTM01 == 0.  This loop defines the initial standing wave on the first tracking step.

        else  // advance the phases smoothly in time and add up the resulting fields again.
        {
        	/*
        FieldFromShort = cos(0.);  // no bounces yet.
        FieldFromPolarizer = cos(0.);  // no bounces yet.
        for (int i=0; i<nbounces; i++)
          {
          phi_shortTM01[i] += 2.*PI*fprime_short*dt;
          phi_polarizerTM01[i] += 2.*PI*fprime_polarizer*dt;
          FieldFromShort += cos(phi_shortTM01[i]);
          FieldFromPolarizer += cos(phi_polarizerTM01[i]);

          }
          */
        }

//    	printf("fieldfromshort is %g\n", FieldFromShort);
//    	printf("fieldfrompolarizer is %g\n", FieldFromPolarizer);
//    	getchar();

  TM01FieldAfterBounces = FieldFromShort + FieldFromPolarizer;

//printf("x y z is %f %f %f\n", aFinalParticle.GetPosition().GetX(), aFinalParticle.GetPosition().GetY(), zposition);
//printf("lambda_short is %g and lambda_polarizer is %g\n", lambda_short, lambda_polarizer);
//printf("TM01FieldAfterBounces is %g\n", TM01FieldAfterBounces); //getchar();

//        printf("TM01FieldAfterBounces is %f\n", TM01FieldAfterBounces);

    	return TM01FieldAfterBounces;
    }



    double CyclotronRadiationExtractor::GetDampingFactor(KSParticle& anInitialParticle, KSParticle& aFinalParticle)
    {
//        double fcyc = aFinalParticle.GetCyclotronFrequency();
//        double GroupVelocity = GetGroupVelocityTE11(aFinalParticle);
//    	double zvelocity = aFinalParticle.GetVelocity().GetZ();
//    	double zposition = aFinalParticle.GetPosition().GetZ();
//        double GammaZ = 1.0/pow(1.0-pow(zvelocity/GetGroupVelocityTE11(aFinalParticle),2.),0.5);

//    	double fprime_short = fcyc*GammaZ*(1.+zvelocity/GroupVelocity);

//    	double phi_short = 2.*PI*2.*(zposition+CENTER_TO_SHORT)/(GroupVelocity/fprime_short);
//        double TE11FieldFromShort = cos(0.) + cos(phi_short); // resonant enhancement.
        double TE11FieldFromShort = GetTE11FieldAfterOneBounce(anInitialParticle, aFinalParticle);
        double TM01FieldAfterBounces = GetTM01FieldAfterBounces(anInitialParticle, aFinalParticle);
        double CouplingFactorTE11 = GetCouplingFactorTE11(aFinalParticle);
        double CouplingFactorTM01 = GetCouplingFactorTM01(aFinalParticle);


        double DampingFactorTE11 = CouplingFactorTE11*(1. - TE11FieldFromShort*TE11FieldFromShort);  // can be > 0 or < 0.
        double DampingFactorTM01 = CouplingFactorTM01*(1. - TM01FieldAfterBounces*TM01FieldAfterBounces);  // can be > 0 or < 0.

//        double DampingFactor = DampingFactorTM01 + DampingFactorTE11;
        double DampingFactor = DampingFactorTE11;



/*
if (fabs(DampingFactor)>0.)
{
//        printf("x, y, z position is %f %f %f\n", aFinalParticle.GetPosition().GetX(), aFinalParticle.GetPosition().GetY(), zposition);
        printf("couplingTM01 is %f\n", CouplingFactorTM01);
        printf("damping factor total = %f, DampingFactorTE11 is %g and DampingFactorTM01 is %g\n", DampingFactor, DampingFactorTE11, DampingFactorTM01);
        getchar();
}
*/

    	return DampingFactor;
    }



    bool CyclotronRadiationExtractor::ExecutePostStepModifcation( KSParticle& anInitialParticle, KSParticle& aFinalParticle, KSParticleQueue& aQueue )
    {


//    	printf("pre step kinetic energy - 4.84338e-15 is %g\n", anInitialParticle.GetKineticEnergy()- 4.84338e-15); //getchar();
//    	printf("post step kinetic energy - 4.84338e-15 is %g\n", aFinalParticle.GetKineticEnergy()- 4.84338e-15); //getchar();


// adjust power with reflections.
    	double DeltaE = GetDampingFactor(anInitialParticle, aFinalParticle)*(aFinalParticle.GetKineticEnergy() - anInitialParticle.GetKineticEnergy());
//    	printf("poststep says DeltaE is %g\n", DeltaE);
    	aFinalParticle.SetKineticEnergy((aFinalParticle.GetKineticEnergy() - DeltaE));

//    	printf("z is %f and DeltaE is %g and post fix kinetic energy is %g\n", aFinalParticle.GetPosition().Z(), DeltaE, aFinalParticle.GetKineticEnergy() - 4.84338e-15); getchar();


    	t_poststep = aFinalParticle.GetTime();

        if (t_poststep - t_old > 5.e-10)
            {
        	std::unique_lock< std::mutex >tLock( fMutexDigitizer, std::defer_lock );  // lock access to mutex before writing to globals.
            tLock.lock();
            Z = aFinalParticle.GetPosition().Z();
            X = aFinalParticle.GetPosition().X();
            Y = aFinalParticle.GetPosition().Y();
            de = aFinalParticle.GetKineticEnergy_eV() - anInitialParticle.GetKineticEnergy_eV();
            dt = aFinalParticle.GetTime() - anInitialParticle.GetTime();
            xvelocity = aFinalParticle.GetVelocity().GetX();
            yvelocity = aFinalParticle.GetVelocity().GetY();
            zvelocity = aFinalParticle.GetVelocity().GetZ();

            GammaZ = 1.0/pow(1.0-pow(zvelocity/GetGroupVelocityTE11(aFinalParticle),2.),0.5);  // speed of light is group velocity

//	    	                fcyc = aFinalParticle.GetCyclotronFrequency();  // inconsistent.  do not use.
					                fcyc = 1.125/dt;
            LarmorPower = -de/dt*1.602677e-19;
            tLock.unlock();
            fDigitizerCondition.notify_one();  // notify Locust after writing.

             t_old = t_poststep;

             /*
             printf("de is %g and dt is %g and LarmorPower is %g\n", de, dt, LarmorPower);

          	 printf("Kassiopeia says:  tick has happened; continuous time is %g and zvelocity is %f\n", t_poststep, zvelocity);
          	 printf("Kassiopeia says:  fcyc is %g\n", fcyc);
          	 printf("  initial particle momentum: %g\n", anInitialParticle.GetMomentum().Z());
             printf("Mass is %g\n", anInitialParticle.GetMass());
             printf("k.e. = %g eV\n", anInitialParticle.GetKineticEnergy_eV());
             printf("Lorentz factor is %f\n", anInitialParticle.GetLorentzFactor());
             printf("1/(sqrt(1-v^2/c^2) is %f\n", 1.0/pow(1.0-pow(anInitialParticle.GetSpeed()/2.99792e8,2.),0.5));
             printf("FinalParticle().IsActive() is %d\n", aFinalParticle.IsActive());
             */
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
