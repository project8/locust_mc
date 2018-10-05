/*
 * LMCFieldCalculator.cc
 *
 *  Created on: Oct 4, 2018
 *      Author: pslocum
 */



#include "LMCFieldCalculator.hh"
#include "LMCGlobalsDeclaration.hh"

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



  double FieldCalculator::GetGroupVelocityTM01(KSParticle& aFinalParticle)
  {
    const double SpeedOfLight = LMCConst::C();
    double CutOffFrequency = 2. * LMCConst::Pi() * SpeedOfLight * 2.405 / 2. / LMCConst::Pi() / 0.00502920; // rad/s
    double cyclotronFrequency = aFinalParticle.GetCyclotronFrequency();
    double GroupVelocity = SpeedOfLight * sqrt( 1. - pow(CutOffFrequency/(2. * LMCConst::Pi() * cyclotronFrequency), 2.) );
    return GroupVelocity;
  }

  double FieldCalculator::GetGroupVelocityTE11(KSParticle& aFinalParticle)
  {
    const double SpeedOfLight = LMCConst::C(); // m/s
    double CutOffFrequency = 2. * LMCConst::Pi() * SpeedOfLight * 1.841 / 2. / LMCConst::Pi() / 0.00502920; // rad/s
    double cyclotronFrequency = aFinalParticle.GetCyclotronFrequency();
    double GroupVelocity = SpeedOfLight * sqrt( 1. - pow(CutOffFrequency / (2. * LMCConst::Pi() * cyclotronFrequency), 2.) );
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



  double FieldCalculator::GetCouplingFactorTE11(KSParticle& aFinalParticle)
  {
    double kc = 1.841/0.00502920;
    double x = aFinalParticle.GetPosition().GetX();
    double y = aFinalParticle.GetPosition().GetY();
        
    double r = sqrt( x * x + y * y);
    double coupling = 119116./168.2 * 2./LMCConst::Pi() * 4./(2.*LMCConst::Pi()) / kc/2. * ( (j0(kc*r) - jn(2,kc*r)) +
											     (j0(kc*r) + jn(2, kc*r)) );
    return coupling*coupling;
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

    //printf("field sum at z=%g is %f with zvelocity %g\n", zPosition, FieldFromShort, zvelocity); 
    //getchar();

    return FieldFromShort;  // Phase 1

  }

  double FieldCalculator::GetTE11FieldAfterOneBounce(KSParticle& anInitialParticle, KSParticle& aFinalParticle)
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

    return TE11FieldAfterOneBounce;
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
    double TE01FieldFromShort = GetTE01FieldAfterOneBounce(anInitialParticle, aFinalParticle);
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
