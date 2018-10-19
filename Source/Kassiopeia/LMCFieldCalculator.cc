/*
 * LMCFieldCalculator.cc
 *
 *  Created on: Oct 4, 2018
 *      Author: pslocum, nbuzinsky
 */



#include "LMCFieldCalculator.hh"
#include "LMCGlobalsDeclaration.hh"
#include <algorithm>

using namespace Kassiopeia;
namespace locust
{

    FieldCalculator::FieldCalculator()
    {
    }

    FieldCalculator::FieldCalculator( const FieldCalculator& aOrig ):  
      Kassiopeia::KSComponent(),
      Kassiopeia::KSComponentTemplate< FieldCalculator, Kassiopeia::KSSpaceInteraction >( aOrig )
    {
    }


    FieldCalculator::~FieldCalculator()
    {
    }


    FieldCalculator* FieldCalculator::Clone() const
    {
        return new FieldCalculator( *this );
    }



  int FieldCalculator::FindNode(double tNew) const
  {
    std::deque<locust::Particle>::iterator it;
    it = std::upper_bound( fParticleHistory.begin() , fParticleHistory.end() , tNew, [] (const double &a , const locust::Particle &b) { return a < b.GetTime();} );

    int tNodeIndex = it - fParticleHistory.begin();

    return tNodeIndex;
  }


  double FieldCalculator::GetSpaceTimeInterval(const double &aParticleTime, const double &aReceiverTime, const LMCThreeVector &aParticlePosition, const LMCThreeVector &aReceiverPosition, double GroupVelocityTE01 )
  {
    //    printf("sti says aReceiverTime is %g, aParticleTime is %g\n", aReceiverTime, aParticleTime);
    return aReceiverTime - aParticleTime - (aReceiverPosition - aParticlePosition).Magnitude() / GroupVelocityTE01;
  }


  double GetFieldStepRoot(const locust::Particle aParticle, double aReceiverTime, LMCThreeVector aReceiverPosition, double aSpaceTimeInterval)
  {
    double tRetardedTime = aParticle.GetTime(true);
    //    printf("tr = %g and  sti = %g\n", tRetardedTime, aSpaceTimeInterval);
    return tRetardedTime + aSpaceTimeInterval;
  }






    double FieldCalculator::GetGroupVelocityTM01(KSParticle& aFinalParticle)
    {
        const double SpeedOfLight = LMCConst::C();
        double CutOffFrequency = 2. * LMCConst::Pi() * SpeedOfLight * 2.405 / 2. / LMCConst::Pi() / 0.00502920; // rad/s
        double cyclotronFrequency = aFinalParticle.GetCyclotronFrequency();
        double GroupVelocity = SpeedOfLight * sqrt( 1. - pow(CutOffFrequency/(2. * LMCConst::Pi() * cyclotronFrequency), 2.) );
        return GroupVelocity;
    }


    double FieldCalculator::GetGroupVelocityTE01(KSParticle& aFinalParticle)  // Phase 1
    {
        double SpeedOfLight = LMCConst::C(); // m/s
        double CutOffFrequency = SpeedOfLight * LMCConst::Pi() / 10.668e-3; // a in m
        double fcyc = aFinalParticle.GetCyclotronFrequency();
        double GroupVelocity = SpeedOfLight * pow( 1. - pow(CutOffFrequency/(2.*LMCConst::Pi()*fcyc), 2.) , 0.5);
        //        printf("GroupVelocity is %g\n", GroupVelocity); getchar();
        return GroupVelocity;
    }


    double FieldCalculator::GetCouplingFactorTM01(KSParticle& aFinalParticle)
    {
        double kc = 2.405/0.00502920;
        double x = aFinalParticle.GetPosition().GetX();
        double y = aFinalParticle.GetPosition().GetY();
        double r = sqrt(x*x + y*y);
        double coupling =   146876.5/168.2 * 2./LMCConst::Pi() * 4./(2.*LMCConst::Pi()) / kc * j1(kc*r);
        return coupling*coupling;
    }


    double FieldCalculator::GetCouplingFactorTE01(KSParticle& aFinalParticle)  // Phase 1
    {
        double dim1_wr42 = 10.668e-3; // a in m
        double x = aFinalParticle.GetPosition().GetX() + dim1_wr42/2.;
        double coupling = 0.63*sin(LMCConst::Pi()*x/dim1_wr42);  // avg over cyclotron orbit.
        return coupling*coupling;
    }


    double FieldCalculator::GetTE01FieldAfterOneBounce(KSParticle& anInitialParticle, KSParticle& aFinalParticle)
    {
        double fcyc = aFinalParticle.GetCyclotronFrequency();
        double GroupVelocity = GetGroupVelocityTE01(aFinalParticle);   
        double zvelocity = aFinalParticle.GetVelocity().GetZ();
        double zPosition = aFinalParticle.GetPosition().GetZ();
        double GammaZ = 1.0/pow(1.0-pow(zvelocity/GetGroupVelocityTE01(aFinalParticle),2.),0.5);

        double fprime_short = fcyc*GammaZ*(1.+zvelocity/GroupVelocity);
        double phi_shortTE01 = LMCConst::Pi()/2. + 2.*LMCConst::Pi()*(fabs(zPosition) + CENTER_TO_SHORT)/(GroupVelocity/fprime_short);  // phase of reflected field at position of electron.          

        double FieldFromShort = cos(0.) + cos(phi_shortTE01);

	//        printf("field sum at z=%g is %f with zvelocity %g\n", zPosition, FieldFromShort, zvelocity); 
        //getchar();

        return FieldFromShort;  // Phase 1

    }


  double FieldCalculator::GetRetardedTE01FieldAfterOneBounce(KSParticle& aFinalParticle)
  {
   
    locust::Particle tCurrentParticle = fParticleHistory.back();
      
    int CurrentIndex;
    
    const double kassiopeiaTimeStep = fabs(fParticleHistory[0].GetTime() - fParticleHistory[1].GetTime());
    const int historySize = fParticleHistory.size();
    double tReceiverTime = aFinalParticle.GetTime();
    double tRetardedTime = 0.; //Retarded time of particle corresponding to when emission occurs, reaching electron at tReceiverTime
    static double fPreviousRetardedTime = 99.;
    static int fPreviousRetardedIndex = -99;
    double tSpaceTimeInterval=0.;
    double dtRetarded=0;
    double tTolerance=1e-23;
    double rxPosition = -2.*CENTER_TO_SHORT - aFinalParticle.GetPosition().GetZ();  // calculate tRetardedTime at grid point on other side of short by -CENTER_TO_SHORT - (z_now - (-CENTER_TO_SHORT))
    LMCThreeVector testPoint(0.,0.,rxPosition);  // electron is at z-now.

    if (tReceiverTime <= fDigitizerTimeStep + 2.*kassiopeiaTimeStep)  // new event starting.  reset.
      {
	fPreviousRetardedTime = 99.;
        fPreviousRetardedIndex = -99;
      }
      

    if (fParticleHistory.front().GetTime()<=3.*kassiopeiaTimeStep)
      {
	fParticleHistory.front().Interpolate(0);
	if(GetSpaceTimeInterval(fParticleHistory.front().GetTime(true), tReceiverTime , fParticleHistory.front().GetPosition(true), testPoint, GetGroupVelocityTE01(aFinalParticle) ) < 0 )
	  {
	    //	    printf("Skipping! out of Bounds!: tReceiverTime=%e\n",tReceiverTime); 
	    //continue;
	  }
      }


    if(fPreviousRetardedIndex == -99)
      {
	CurrentIndex=FindNode(tReceiverTime);
	tCurrentParticle = fParticleHistory[CurrentIndex];
	tRetardedTime = tReceiverTime - (tCurrentParticle.GetPosition() - testPoint ).Magnitude() /  GetGroupVelocityTE01(aFinalParticle);
      }
    else
      {
	CurrentIndex = fPreviousRetardedIndex;
	tRetardedTime = fPreviousRetardedTime + kassiopeiaTimeStep;
      }

    CurrentIndex = FindNode(tRetardedTime);
    CurrentIndex = std::min(std::max(CurrentIndex,0) , historySize - 1);

    tCurrentParticle = fParticleHistory[CurrentIndex];
    tCurrentParticle.Interpolate(tRetardedTime);
    tSpaceTimeInterval = GetSpaceTimeInterval(tCurrentParticle.GetTime(true), tReceiverTime, tCurrentParticle.GetPosition(true), testPoint, GetGroupVelocityTE01(aFinalParticle));

    double tOldSpaceTimeInterval=99.;

    //Converge to root
    for(int j=0;j<25;++j)
      {
	tRetardedTime = GetFieldStepRoot(tCurrentParticle, tReceiverTime, testPoint, tSpaceTimeInterval);
	tCurrentParticle.Interpolate(tRetardedTime);

	//Change the kassiopeia step we expand around if the interpolation time displacement is too large
	if(fabs(tCurrentParticle.GetTime(true) - tCurrentParticle.GetTime(false)) > kassiopeiaTimeStep)
	  {
	    CurrentIndex=FindNode(tRetardedTime);
	    tCurrentParticle=fParticleHistory[CurrentIndex];
	    tCurrentParticle.Interpolate(tRetardedTime);
	  }

	tSpaceTimeInterval = GetSpaceTimeInterval(tCurrentParticle.GetTime(true), tReceiverTime, tCurrentParticle.GetPosition(true), testPoint, GetGroupVelocityTE01(aFinalParticle));
	tOldSpaceTimeInterval = tSpaceTimeInterval;
      }


    fPreviousRetardedIndex = CurrentIndex;
    fPreviousRetardedTime = tRetardedTime;


      

    double fcyc = tCurrentParticle.GetCyclotronFrequency()/2./LMCConst::Pi();
    double GroupVelocity = GetGroupVelocityTE01(aFinalParticle);
    double zvelocity = tCurrentParticle.GetVelocity().Z();  // v retarded particle.
    double zPosition = aFinalParticle.GetPosition().GetZ();  // z now.
    double GammaZ = 1.0/pow(1.0-pow(zvelocity/GroupVelocity,2.),0.5);

    double fprime_short = fcyc*GammaZ*(1.+zvelocity/GroupVelocity);
    double phi_shortTE01 = LMCConst::Pi()/2. + 2.*LMCConst::Pi()*(fabs(zPosition) + CENTER_TO_SHORT)/(GroupVelocity/fprime_short);  // phase of reflected field at position of electron.  force mode symmetry (critical) with fabs(z).
    double FieldFromShort = cos(0.) + cos(phi_shortTE01);


    //    printf("currentIndex is %d and zvelocity is %g and fPreviousRetardedIndex is %d and fPreviousRetardedTime is %g\n", CurrentIndex, zvelocity, fPreviousRetardedIndex, fPreviousRetardedTime); getchar();
    //printf("field sum at z=%g is %f with zvelocity %g and CurrentIndex %d and fcyc is %g\n", zPosition, FieldFromShort, zvelocity, CurrentIndex, fcyc);
    //    printf("particle time is %g and t_old is %g\n", aFinalParticle.GetTime(), t_old);
    //    getchar();                                                           
   
return FieldFromShort;

}



    double FieldCalculator::GetTM01FieldWithTerminator(KSParticle& anInitialParticle, KSParticle& aFinalParticle)
    {
        double tCyclotronFrequency = aFinalParticle.GetCyclotronFrequency();
        double GroupVelocity = GetGroupVelocityTM01(aFinalParticle);
        double tVelocityZ = aFinalParticle.GetVelocity().GetZ();
        double tPositionZ = aFinalParticle.GetPosition().GetZ();
        double GammaZ = 1.0/sqrt(1.0-pow(tVelocityZ / GetGroupVelocityTM01(aFinalParticle),2.) );
        double fprime_polarizer = tCyclotronFrequency*GammaZ*(1.-tVelocityZ/GroupVelocity);

        double phi_polarizerTM01 = 2.*LMCConst::Pi()*(2.*(CENTER_TO_ANTENNA-fabs(tPositionZ)))/(GroupVelocity/fprime_polarizer);
        double TM01FieldWithTerminator = cos(0.) + cos(phi_polarizerTM01);
        //printf("TM01FieldWithTerminator is %f\n", TM01FieldWithTerminator);
        //getchar();
        return TM01FieldWithTerminator;


    }


    double FieldCalculator::GetDampingFactorPhase1(KSParticle& anInitialParticle, KSParticle& aFinalParticle)
    {
      double TE01FieldFromShort = 0.;
      if (fParticleHistory.size()&&(0==0)) // if fParticleHistory has some entries.
	{
        TE01FieldFromShort = GetRetardedTE01FieldAfterOneBounce(aFinalParticle);
        }
      else  // if fParticleHistory has not been filled at all yet.
        {
	TE01FieldFromShort = GetTE01FieldAfterOneBounce(anInitialParticle, aFinalParticle);
	}

      //      if (fParticleHistory.size()) {printf("FieldFromShort is %g\n", TE01FieldFromShort); getchar();}


        double A10squ = GetCouplingFactorTE01(aFinalParticle);
        double DampingFactorTE01 = 1. - A10squ + A10squ*TE01FieldFromShort*TE01FieldFromShort;  // = P'/P

        return DampingFactorTE01;
    }




    double FieldCalculator::GetDampingFactorPhase2(KSParticle& anInitialParticle, KSParticle& aFinalParticle)
    {
        double TM01FieldWithTerminator = GetTM01FieldWithTerminator(anInitialParticle, aFinalParticle);
        double A01squ = GetCouplingFactorTM01(aFinalParticle);
        double DampingFactorTM01 = 1. - A01squ + A01squ*TM01FieldWithTerminator*TM01FieldWithTerminator;  // = P'/P
        double DampingFactor = DampingFactorTM01;

        return DampingFactor;
    }








    void FieldCalculator::CalculateInteraction(
            const Kassiopeia::KSTrajectory& aTrajectory,
            const Kassiopeia::KSParticle& aTrajectoryInitialParticle,
            const Kassiopeia::KSParticle& aTrajectoryFinalParticle,
            const KThreeVector& aTrajectoryCenter,
            const double& aTrajectoryRadius,
            const double& aTrajectoryTimeStep,
            Kassiopeia::KSParticle& anInteractionParticle,
            double& anInteractionStep, bool& anInteractionFlag
            )
    {
    }

    void FieldCalculator::ExecuteInteraction(
            const Kassiopeia::KSParticle& anInitialParticle,
            Kassiopeia::KSParticle& aFinalParticle,
            Kassiopeia::KSParticleQueue& aSecondaries
            ) const
    {
    }




}  /* namespace locust */
