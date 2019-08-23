
#include "LMCFieldCalculator.hh"
#include "KSParticleFactory.h"
#include <algorithm>

#include "KSInteractionsMessage.h"
#include <limits>

using std::numeric_limits;

using namespace Kassiopeia;
namespace locust
{

    FieldCalculator::FieldCalculator() :
            fSpaceInteractions( 128 ),
            fSpaceInteraction( NULL ),
            fStep( NULL ),
            fTerminatorParticle( NULL ),
            fTrajectoryParticle( NULL ),
            fInteractionParticle( NULL ),
            fFinalParticle( NULL ),
            fParticleQueue( NULL ),
            fTrajectory( NULL ),
            fInterface( KLInterfaceBootstrapper::get_instance()->GetInterface() )
    {
    }
    FieldCalculator::FieldCalculator( const FieldCalculator& aCopy ) :
            KSComponent(),
            fSpaceInteractions( aCopy.fSpaceInteractions ),
            fSpaceInteraction( aCopy.fSpaceInteraction ),
            fStep( aCopy.fStep ),
            fTerminatorParticle( aCopy.fTerminatorParticle ),
            fTrajectoryParticle( aCopy.fTrajectoryParticle ),
            fInteractionParticle( aCopy.fInteractionParticle ),
            fFinalParticle( aCopy.fFinalParticle ),
            fParticleQueue( aCopy.fParticleQueue ),
            fTrajectory( aCopy.fTrajectory ),
            fInterface( aCopy.fInterface )
    {
    }
    FieldCalculator* FieldCalculator::Clone() const
    {
        return new FieldCalculator( *this );
    }
    FieldCalculator::~FieldCalculator()
    {
    }



    double FieldCalculator::GetGroupVelocityTM01(KSParticle& aFinalParticle)
    {
        const double SpeedOfLight = LMCConst::C();
        double CutOffFrequency = 2. * LMCConst::Pi() * SpeedOfLight * 2.405 / 2. / LMCConst::Pi() / 0.00502920; // rad/s
        double cyclotronFrequency = aFinalParticle.GetCyclotronFrequency();
        double GroupVelocity = SpeedOfLight * sqrt( 1. - pow(CutOffFrequency/(2. * LMCConst::Pi() * cyclotronFrequency), 2.) );
        return GroupVelocity;
    }


    double FieldCalculator::GetGroupVelocityTE10(KSParticle& aFinalParticle)  // Phase 1
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

    // power fraction. 48.1 is numerical normalization
   // of J \cdot E after time averaging as in Collin IEEE paper.
   // TM power reduction of 0.00136 is included in normalization as in Collin paper.
   // coupling*coupling is the power fraction plotted in the Locust paper.

        double coupling =   705.7 * 2./LMCConst::Pi() * 4./(2.*LMCConst::Pi()) / kc *
        		j1(kc*r);
      
        return coupling*coupling;
    }


    double FieldCalculator::GetCouplingFactorTE10(KSParticle& aFinalParticle)  // Phase 1
    {
        double dim1_wr42 = 10.668e-3; // a in m
        double x = aFinalParticle.GetPosition().GetX() + dim1_wr42/2.;
        double coupling = 0.63*sin(LMCConst::Pi()*x/dim1_wr42);  // avg over cyclotron orbit.
        return coupling*coupling;
    }


    double FieldCalculator::GetTE10FieldAfterOneBounce(KSParticle& aFinalParticle)
    {
        double fcyc = aFinalParticle.GetCyclotronFrequency();
        double GroupVelocity = GetGroupVelocityTE10(aFinalParticle);
        double zvelocity = aFinalParticle.GetVelocity().GetZ();
        double zPosition = aFinalParticle.GetPosition().GetZ();
        double GammaZ = 1.0/pow(1.0-pow(zvelocity/GetGroupVelocityTE10(aFinalParticle),2.),0.5);

        double fprime_short = fcyc*GammaZ*(1.+zvelocity/GroupVelocity);
        double phi_shortTE10 = LMCConst::Pi()/2. + 2.*LMCConst::Pi()*(fabs(zPosition) + fInterface->fCENTER_TO_SHORT)/(GroupVelocity/fprime_short);  // phase of reflected field at position of electron.

        double FieldFromShort = cos(0.) + cos(phi_shortTE10);

	//        printf("field sum at z=%g is %f with zvelocity %g\n", zPosition, FieldFromShort, zvelocity); 
        //getchar();

        return FieldFromShort;  // Phase 1

    }



    double FieldCalculator::GetTM01FieldWithTerminator(KSParticle& aFinalParticle)
    {
        double tCyclotronFrequency = aFinalParticle.GetCyclotronFrequency();
        double GroupVelocity = GetGroupVelocityTM01(aFinalParticle);
        double tVelocityZ = aFinalParticle.GetVelocity().GetZ();
        double tPositionZ = aFinalParticle.GetPosition().GetZ();
        double GammaZ = 1.0/sqrt(1.0-pow(tVelocityZ / GetGroupVelocityTM01(aFinalParticle),2.) );
        double fprime_polarizer = tCyclotronFrequency*GammaZ*(1.-tVelocityZ/GroupVelocity);

        double phi_polarizerTM01 = 2.*LMCConst::Pi()*(2.*(fInterface->fCENTER_TO_ANTENNA-fabs(tPositionZ)))/(GroupVelocity/fprime_polarizer);
        double TM01FieldWithTerminator = cos(0.) + cos(phi_polarizerTM01);
        //printf("TM01FieldWithTerminator is %f\n", TM01FieldWithTerminator);
        //getchar();
        return TM01FieldWithTerminator;


    }


    double FieldCalculator::GetDampingFactorPhase1(KSParticle& aFinalParticle)
    {
      double TE10FieldFromShort = 0.;
	  TE10FieldFromShort = GetTE10FieldAfterOneBounce(aFinalParticle);
      double A10squ = GetCouplingFactorTE10(aFinalParticle);
      double DampingFactorTE10 = 1. - A10squ + A10squ*TE10FieldFromShort*TE10FieldFromShort;  // = P'/P
      return DampingFactorTE10;
    }



    double FieldCalculator::GetDampingFactorPhase2(KSParticle& aFinalParticle)
    {
      double TM01FieldWithTerminator = 0.;
      TM01FieldWithTerminator = GetTM01FieldWithTerminator(aFinalParticle);
      double A01squ = GetCouplingFactorTM01(aFinalParticle);
      double DampingFactorTM01 = 1. - A01squ + A01squ*TM01FieldWithTerminator*TM01FieldWithTerminator;  // = P'/P
      return DampingFactorTM01;
    }






    void FieldCalculator::CalculateInteraction( const KSTrajectory& aTrajectory, const KSParticle& aTrajectoryInitialParticle, const KSParticle& aTrajectoryFinalParticle, const KThreeVector& aTrajectoryCenter, const double& aTrajectoryRadius, const double& aTrajectoryStep, KSParticle& anInteractionParticle, double& anInteractionStep, bool& anInteractionFlag )
    {
        KSParticle tInteractionParticle;
        double tInteractionStep;
        bool tInteractionFlag;

        anInteractionParticle = aTrajectoryFinalParticle;
        anInteractionStep = aTrajectoryStep;
        anInteractionFlag = false;
        for( int tIndex = 0; tIndex < fSpaceInteractions.End(); tIndex++ )
        {
            fSpaceInteractions.ElementAt( tIndex )->CalculateInteraction( aTrajectory, aTrajectoryInitialParticle, aTrajectoryFinalParticle, aTrajectoryCenter, aTrajectoryRadius, aTrajectoryStep, tInteractionParticle, tInteractionStep, tInteractionFlag );
            if( tInteractionFlag == true )
            {
                anInteractionFlag = true;
                if( tInteractionStep < anInteractionStep )
                {
                    anInteractionParticle = tInteractionParticle;
                    fSpaceInteraction = fSpaceInteractions.ElementAt( tIndex );
                }
            }
        }
        return;
    }
    void FieldCalculator::ExecuteInteraction( const KSParticle& anInteractionParticle, KSParticle& aFinalParticle, KSParticleQueue& aSecondaries ) const
    {
        if( fSpaceInteraction != NULL )
        {
            fSpaceInteraction->ExecuteInteraction( anInteractionParticle, aFinalParticle, aSecondaries );
        }
        else
        {
            aFinalParticle = anInteractionParticle;
        }
        return;
    }

    void FieldCalculator::AddSpaceInteraction( KSSpaceInteraction* aSpaceInteraction )
    {
        if( fSpaceInteractions.AddElement( aSpaceInteraction ) == -1 )
        {
            intmsg( eError ) << "<" << GetName() << "> could not add space interaction <" << aSpaceInteraction->GetName() << ">" << eom;
            return;
        }
        intmsg_debug( "<" << GetName() << "> adding space interaction <" << aSpaceInteraction->GetName() << ">" << eom )
        return;
    }
    void FieldCalculator::RemoveSpaceInteraction( KSSpaceInteraction* aSpaceInteraction )
    {
        if( fSpaceInteractions.RemoveElement( aSpaceInteraction ) == -1 )
        {
            intmsg( eError ) << "<" << GetName() << "> could not remove space interaction <" << aSpaceInteraction->GetName() << ">" << eom;
            return;
        }
        intmsg_debug( "<" << GetName() << "> removing space interaction <" << aSpaceInteraction->GetName() << ">" << eom )
        return;
    }

    void FieldCalculator::SetStep( KSStep* aStep )
    {
        fStep = aStep;
        fTerminatorParticle = &(aStep->TerminatorParticle());
        fTrajectoryParticle = &(aStep->TrajectoryParticle());
        fInteractionParticle = &(aStep->InteractionParticle());
        fFinalParticle = &(aStep->FinalParticle());
        fParticleQueue = &(aStep->ParticleQueue());
        return;
    }
    void FieldCalculator::SetTrajectory( KSTrajectory* aTrajectory )
    {
        fTrajectory = aTrajectory;
        return;
    }

    void FieldCalculator::CalculateInteraction()
    {
        *fInteractionParticle = *fTrajectoryParticle;

        if( fSpaceInteractions.End() == 0 )
        {
            intmsg_debug( "space interaction calculation:" << eom )
            intmsg_debug( "  no space interactions active" << eom )
            intmsg_debug( "  interaction name: <" << fStep->GetSpaceInteractionName() << ">" << eom )
            intmsg_debug( "  interaction step: <" << fStep->GetSpaceInteractionStep() << ">" << eom )
            intmsg_debug( "  interaction flag: <" << fStep->GetSpaceInteractionFlag() << ">" << eom )

            intmsg_debug( "space interaction calculation interaction particle state: " << eom )
            intmsg_debug( "  final particle space: <" << (fInteractionParticle->GetCurrentSpace() ? fInteractionParticle->GetCurrentSpace()->GetName() : "" ) << ">" << eom )
            intmsg_debug( "  final particle surface: <" << (fInteractionParticle->GetCurrentSurface() ? fInteractionParticle->GetCurrentSurface()->GetName() : "" ) << ">" << eom )
            intmsg_debug( "  final particle time: <" << fInteractionParticle->GetTime() << ">" << eom )
            intmsg_debug( "  final particle length: <" << fInteractionParticle->GetLength() << ">" << eom )
            intmsg_debug( "  final particle position: <" << fInteractionParticle->GetPosition().X() << ", " << fInteractionParticle->GetPosition().Y() << ", " << fInteractionParticle->GetPosition().Z() << ">" << eom )
            intmsg_debug( "  final particle momentum: <" << fInteractionParticle->GetMomentum().X() << ", " << fInteractionParticle->GetMomentum().Y() << ", " << fInteractionParticle->GetMomentum().Z() << ">" << eom )
            intmsg_debug( "  final particle kinetic energy: <" << fInteractionParticle->GetKineticEnergy_eV() << ">" << eom )
            intmsg_debug( "  final particle electric field: <" << fInteractionParticle->GetElectricField().X() << "," << fInteractionParticle->GetElectricField().Y() << "," << fInteractionParticle->GetElectricField().Z() << ">" << eom )
            intmsg_debug( "  final particle magnetic field: <" << fInteractionParticle->GetMagneticField().X() << "," << fInteractionParticle->GetMagneticField().Y() << "," << fInteractionParticle->GetMagneticField().Z() << ">" << eom )
            intmsg_debug( "  final particle angle to magnetic field: <" << fInteractionParticle->GetPolarAngleToB() << ">" << eom )

            return;
        }

        CalculateInteraction( *fTrajectory, *fTerminatorParticle, *fTrajectoryParticle, fStep->TrajectoryCenter(), fStep->TrajectoryRadius(), fStep->TrajectoryStep(), *fInteractionParticle, fStep->SpaceInteractionStep(), fStep->SpaceInteractionFlag() );

        if( fStep->SpaceInteractionFlag() == true )
        {
            intmsg_debug( "space interaction calculation:" << eom )
            intmsg_debug( "  space interaction may occur" << eom )
        }
        else
        {
            intmsg_debug( "space interaction calculation:" << eom )
            intmsg_debug( "  space interaction will not occur" << eom )
        }

        intmsg_debug( "space interaction calculation interaction particle state: " << eom )
        intmsg_debug( "  interaction particle space: <" << (fInteractionParticle->GetCurrentSpace() ? fInteractionParticle->GetCurrentSpace()->GetName() : "" ) << ">" << eom )
        intmsg_debug( "  interaction particle surface: <" << (fInteractionParticle->GetCurrentSurface() ? fInteractionParticle->GetCurrentSurface()->GetName() : "" ) << ">" << eom )
        intmsg_debug( "  interaction particle time: <" << fInteractionParticle->GetTime() << ">" << eom )
        intmsg_debug( "  interaction particle length: <" << fInteractionParticle->GetLength() << ">" << eom )
        intmsg_debug( "  interaction particle position: <" << fInteractionParticle->GetPosition().X() << ", " << fInteractionParticle->GetPosition().Y() << ", " << fInteractionParticle->GetPosition().Z() << ">" << eom )
        intmsg_debug( "  interaction particle momentum: <" << fInteractionParticle->GetMomentum().X() << ", " << fInteractionParticle->GetMomentum().Y() << ", " << fInteractionParticle->GetMomentum().Z() << ">" << eom )
        intmsg_debug( "  interaction particle kinetic energy: <" << fInteractionParticle->GetKineticEnergy_eV() << ">" << eom )
        intmsg_debug( "  interaction particle electric field: <" << fInteractionParticle->GetElectricField().X() << "," << fInteractionParticle->GetElectricField().Y() << "," << fInteractionParticle->GetElectricField().Z() << ">" << eom )
        intmsg_debug( "  interaction particle magnetic field: <" << fInteractionParticle->GetMagneticField().X() << "," << fInteractionParticle->GetMagneticField().Y() << "," << fInteractionParticle->GetMagneticField().Z() << ">" << eom )
        intmsg_debug( "  interaction particle angle to magnetic field: <" << fInteractionParticle->GetPolarAngleToB() << ">" << eom );

        return;
    }

    void FieldCalculator::ExecuteInteraction()
    {
        ExecuteInteraction( *fInteractionParticle, *fFinalParticle, *fParticleQueue );
        fFinalParticle->ReleaseLabel( fStep->SpaceInteractionName() );

        fStep->ContinuousTime() = fInteractionParticle->GetTime() - fTerminatorParticle->GetTime();
        fStep->ContinuousLength() = fInteractionParticle->GetLength() - fTerminatorParticle->GetLength();
        fStep->ContinuousEnergyChange() = fInteractionParticle->GetKineticEnergy_eV() - fTerminatorParticle->GetKineticEnergy_eV();
        fStep->ContinuousMomentumChange() = (fInteractionParticle->GetMomentum() - fTerminatorParticle->GetMomentum()).Magnitude();
        fStep->DiscreteSecondaries() = fParticleQueue->size();
        fStep->DiscreteEnergyChange() = fFinalParticle->GetKineticEnergy_eV() - fInteractionParticle->GetKineticEnergy_eV();
        fStep->DiscreteMomentumChange() = (fFinalParticle->GetMomentum() - fInteractionParticle->GetMomentum()).Magnitude();

        intmsg_debug( "space interaction execution:" << eom )
        intmsg_debug( "  space interaction name: <" << fStep->SpaceInteractionName() << ">" << eom )
        intmsg_debug( "  step continuous time: <" << fStep->ContinuousTime() << ">" << eom )
        intmsg_debug( "  step continuous length: <" << fStep->ContinuousLength() << ">" << eom )
        intmsg_debug( "  step continuous energy change: <" << fStep->ContinuousEnergyChange() << ">" << eom )
        intmsg_debug( "  step continuous momentum change: <" << fStep->ContinuousMomentumChange() << ">" << eom )
        intmsg_debug( "  step discrete secondaries: <" << fStep->DiscreteSecondaries() << ">" << eom )
        intmsg_debug( "  step discrete energy change: <" << fStep->DiscreteEnergyChange() << ">" << eom )
        intmsg_debug( "  step discrete momentum change: <" << fStep->DiscreteMomentumChange() << ">" << eom );

        intmsg_debug( "space interaction execution final particle state: " << eom )
        intmsg_debug( "  final particle space: <" << (fFinalParticle->GetCurrentSpace() ? fFinalParticle->GetCurrentSpace()->GetName() : "" ) << ">" << eom )
        intmsg_debug( "  final particle surface: <" << (fFinalParticle->GetCurrentSurface() ? fFinalParticle->GetCurrentSurface()->GetName() : "" ) << ">" << eom )
        intmsg_debug( "  final particle time: <" << fFinalParticle->GetTime() << ">" << eom )
        intmsg_debug( "  final particle length: <" << fFinalParticle->GetLength() << ">" << eom )
        intmsg_debug( "  final particle position: <" << fFinalParticle->GetPosition().X() << ", " << fFinalParticle->GetPosition().Y() << ", " << fFinalParticle->GetPosition().Z() << ">" << eom )
        intmsg_debug( "  final particle momentum: <" << fFinalParticle->GetMomentum().X() << ", " << fFinalParticle->GetMomentum().Y() << ", " << fFinalParticle->GetMomentum().Z() << ">" << eom )
        intmsg_debug( "  final particle kinetic energy: <" << fFinalParticle->GetKineticEnergy_eV() << ">" << eom )
        intmsg_debug( "  final particle electric field: <" << fFinalParticle->GetElectricField().X() << "," << fFinalParticle->GetElectricField().Y() << "," << fFinalParticle->GetElectricField().Z() << ">" << eom )
        intmsg_debug( "  final particle magnetic field: <" << fFinalParticle->GetMagneticField().X() << "," << fFinalParticle->GetMagneticField().Y() << "," << fFinalParticle->GetMagneticField().Z() << ">" << eom )
        intmsg_debug( "  final particle angle to magnetic field: <" << fFinalParticle->GetPolarAngleToB() << ">" << eom );

        return;
    }

    void FieldCalculator::PushUpdateComponent()
    {
        for( int tIndex = 0; tIndex < fSpaceInteractions.End(); tIndex++ )
        {
            fSpaceInteractions.ElementAt( tIndex )->PushUpdate();
        }
    }
    void FieldCalculator::PushDeupdateComponent()
    {
        for( int tIndex = 0; tIndex < fSpaceInteractions.End(); tIndex++ )
        {
            fSpaceInteractions.ElementAt( tIndex )->PushDeupdate();
        }
    }

    STATICINT sFieldCalculatorDict = KSDictionary< FieldCalculator >::AddCommand( &FieldCalculator::AddSpaceInteraction, &FieldCalculator::RemoveSpaceInteraction, "add_space_interaction", "remove_space_interaction" );

}
