/*
 * LMCParticleSlim.cc
 *
 *  Created on: Jan 13, 2017
 *      Author: nbuzinsky
 */

#include "LMCParticleSlim.hh"


namespace locust
{

    ParticleSlim::ParticleSlim() :
            fTime( -99. ), 
            fTimeDisplacement( -99.), 
            fTimeStep( -99.), 
            fPosition( -99., -99., -99. ),
            fVelocity( -99., -99., -99. ),
            fVelocityParallel( -99.),
            fNewPosition( -99., -99., -99. ),
            fNewVelocity( -99., -99., -99. ),
            fGuidingCenterPosition( -99., -99., -99. ),
            fMagneticField( -99., -99., -99. ),
            fAlpha( -99., -99., -99. ),
            fBeta( -99., -99., -99. ),
            fMass( -99. ),
            fCharge( -99. ),
            fGamma( -99. ),
            fCyclotronFrequency( -99. ),
            fCyclotronRadius( -99. ),
            fSplineC( -99., -99., -99. ),
            fSplineD( -99., -99., -99. ),
            fReceiverTime( -99. ),
            fReceiverPosition( -99., -99., -99. ),
            fReceiverVector( -99., -99., -99. ),
            fReceiverDistance( -99. ),
            fReceiverDir( -99., -99., -99. )
    {
    }

    ParticleSlim::~ParticleSlim()
    {

    }

    void ParticleSlim::SetTime(double aTime)
    {
        fTime = aTime;
        fTimeDisplacement = 0.;
        return;
    }
    double ParticleSlim::GetTime()
    {
        return fTime;
    }
    double ParticleSlim::GetTimeDisplacement()
    {
        return fTimeDisplacement;
    }

    void ParticleSlim::SetPosition(double X, double Y, double Z)
    {
        fPosition.SetComponents(X,Y,Z);
        fNewPosition.SetComponents(X,Y,Z);
        return;
    }

    KGeoBag::KThreeVector ParticleSlim::GetPosition()
    {
        return fPosition;
    }

    void ParticleSlim::SetVelocityVector(double Vx, double Vy, double Vz)
    {
        fVelocity.SetComponents(Vx,Vy,Vz);
        fNewVelocity.SetComponents(Vx,Vy,Vz);
        return;
    }

    KGeoBag::KThreeVector ParticleSlim::GetVelocity()
    {
        return fVelocity;
    }


    void ParticleSlim::SetMagneticFieldVector(double Bx, double By, double Bz)
    {
        fMagneticField.SetComponents(Bx,By,Bz);
        return;
    }
    void ParticleSlim::SetCyclotronFrequency(double fcyc)
    {
        fCyclotronFrequency=fcyc;
        return;
    }

    void ParticleSlim::SetMass(double aMass)
    {
        fMass = aMass;
        return;
    }

    void ParticleSlim::SetCharge(double aCharge)
    {
        fCharge = aCharge;
        return;
    }

    void ParticleSlim::SetReceiverPosition(double X, double Y, double Z)
    {
        fReceiverPosition.SetComponents(X,Y,Z);
        fReceiverVector=fReceiverPosition-fNewPosition;
        fReceiverDistance=fReceiverVector.Magnitude();
        fReceiverDir=fReceiverVector.Unit();
        return;
    }

    void ParticleSlim::SetReceiverTime(double aTime)
    {
        fReceiverTime=aTime;
        return;
    }

    void ParticleSlim::SetKinematicProperties()
    {
        double beta=fVelocity.Magnitude()/KConst::C();
        fGamma=1./sqrt((1.-beta)*(1.+beta));

        //fCyclotronFrequency=fabs(fCharge) * fMagneticField.Magnitude() / ( fMass * fGamma);
        //printf("%f8.4\n",fCyclotronFrequency);
        fVelocityParallel=fVelocity.Dot(fMagneticField.Unit());
        KGeoBag::KThreeVector vPerp=fVelocity-fVelocityParallel*fMagneticField.Unit();
        fBeta=vPerp.Unit();

        fCyclotronRadius = fGamma * fMass * vPerp.Magnitude() / (fabs(fCharge) * fMagneticField.Magnitude());

        KGeoBag::KThreeVector tCrossRadial=(fVelocity.Cross(fMagneticField)).Unit();

        if(fCharge<0)
        {
            tCrossRadial=-1.*tCrossRadial;
        }

        fGuidingCenterPosition=fPosition+fCyclotronRadius*tCrossRadial;

        fAlpha=fPosition-fGuidingCenterPosition;

        //fBeta=fAlpha.Cross(fMagneticField);
        fAlpha=fAlpha.Unit();
        //fBeta=vPerp.Unit();

        fSplineC.SetComponents(0.,0.,0.);
        fSplineD.SetComponents(0.,0.,0.);

        return;
    }

    void ParticleSlim::Interpolate(double tNew)
    {
        fTimeDisplacement=tNew-fTime; 
        
        double CyclotronFrequency=fCyclotronFrequency;
        if(fCharge<0)CyclotronFrequency*=-1.;

        fNewPosition=fGuidingCenterPosition + fVelocityParallel * fMagneticField.Unit() * fTimeDisplacement + fCyclotronRadius * ( cos(CyclotronFrequency*fTimeDisplacement)*fAlpha + sin(CyclotronFrequency*fTimeDisplacement)*fBeta)+fSplineC*pow(fTimeDisplacement/fTimeStep,3.)+fSplineD*pow(fTimeDisplacement/fTimeStep,4.);

        fNewVelocity=fVelocityParallel*fMagneticField.Unit() + fCyclotronRadius*fCyclotronFrequency*(-sin(CyclotronFrequency*fTimeDisplacement)*fAlpha+cos(CyclotronFrequency*fTimeDisplacement)*fBeta)+3.*fSplineC*pow(fTimeDisplacement/fTimeStep,2.)/fTimeStep+4.*fSplineD*pow(fTimeDisplacement/fTimeStep,3.)/fTimeStep;

        return;
    }

    double ParticleSlim::GetSpaceTimeInterval()
    {
        double s=fReceiverTime-fTime-fTimeDisplacement-fReceiverVector.Magnitude()/KConst::C();
        return s;
    }

    double ParticleSlim::GetSpaceTimeInterval(double aTime)
    {
        double dt=aTime-fTime;

        double CyclotronFrequency=fCyclotronFrequency;
        if(fCharge<0)CyclotronFrequency*=-1.;

        KGeoBag::KThreeVector newPos = fGuidingCenterPosition + fVelocityParallel * fMagneticField.Unit() * dt + fCyclotronRadius * ( cos(CyclotronFrequency*dt)*fAlpha + sin(CyclotronFrequency*dt)*fBeta) + fSplineC *pow(dt/fTimeStep,3.) + fSplineD*pow(dt/fTimeStep,4.);

        KGeoBag::KThreeVector dr=fReceiverPosition-newPos;
        double s=fReceiverTime-aTime-dr.Magnitude()/KConst::C();
        return s;
    }


    double ParticleSlim::GetReceiverDistance()
    {
        return fReceiverDistance;
    }

  //  double ParticleSlim::NewtonStep(double dtRetarded)
  //  {
  //      double tRetarded=fTime+fTimeDisplacement;

  //      double A,B,C;
  //      CalculateQuadraticCoefficients(A,B,C);

  //      double c=KConst::C();

  //      double f_sin=A*dtRetarded*dtRetarded+B*dtRetarded+C+2.*fCyclotronRadius/(c*c)*(fReceiverVector.Dot((cos(fCyclotronFrequency*(dtRetarded+tRetarded))-cos(fCyclotronFrequency*tRetarded))*fAlpha+(sin(fCyclotronFrequency*(dtRetarded+tRetarded))-sin(fCyclotronFrequency*tRetarded))*fBeta)-fCyclotronRadius*(1.-cos(fCyclotronFrequency*dtRetarded)));
  //      double df_sin=2.*A*dtRetarded+B+2.*fCyclotronRadius/(c*c)*fCyclotronFrequency*(fReceiverVector.Dot(-sin(fCyclotronFrequency*(dtRetarded+tRetarded))*fAlpha+cos(fCyclotronFrequency*(dtRetarded+tRetarded))*fBeta)-fCyclotronRadius*sin(fCyclotronFrequency*dtRetarded));

  //      return dtRetarded-f_sin/df_sin;
  //  }

  //  double ParticleSlim::HouseHolderStep(double dtRetarded)
  //  {
  //      double tRetarded=fTime+fTimeDisplacement;

  //      double A,B,C;
  //      CalculateQuadraticCoefficients(A,B,C);

  //      double c=KConst::C();

  //      double f_sin=A*dtRetarded*dtRetarded+B*dtRetarded+C+2.*fCyclotronRadius/(c*c)*(fReceiverVector.Dot((cos(fCyclotronFrequency*(dtRetarded+tRetarded))-cos(fCyclotronFrequency*tRetarded))*fAlpha+(sin(fCyclotronFrequency*(dtRetarded+tRetarded))-sin(fCyclotronFrequency*tRetarded))*fBeta)-fCyclotronRadius*(1.-cos(fCyclotronFrequency*dtRetarded)));
  //      double df_sin=2.*A*dtRetarded+B+2.*fCyclotronRadius/(c*c)*fCyclotronFrequency*(fReceiverVector.Dot(-sin(fCyclotronFrequency*(dtRetarded+tRetarded))*fAlpha+cos(fCyclotronFrequency*(dtRetarded+tRetarded))*fBeta)-fCyclotronRadius*sin(fCyclotronFrequency*dtRetarded));
  //      double d2f_sin=2.*A - 2.*fCyclotronRadius/(c*c)*fCyclotronFrequency*fCyclotronFrequency*(fReceiverVector.Dot(cos(fCyclotronFrequency*(dtRetarded+tRetarded))*fAlpha+sin(fCyclotronFrequency*(dtRetarded+tRetarded))*fBeta)+fCyclotronRadius*cos(fCyclotronFrequency*dtRetarded));

  //      double SqrtArg=1. - 2. * d2f_sin * f_sin / (df_sin * df_sin);
  //      if(SqrtArg < 0)  return NewtonStep(dtRetarded);

  //      return dtRetarded - df_sin / d2f_sin * ( 1. - sqrt(SqrtArg));
  //  }


    double ParticleSlim::CalculatePower(double X, double Y, double Z)
    {
        double c=KConst::C();

        KGeoBag::KThreeVector SurfaceNormal(X,Y,Z);

        //Lorentz Equation
        KGeoBag::KThreeVector BetaDot=fCharge/(fMass * c) *fNewVelocity.Cross(fMagneticField);

        //Lienard-Wiechert Equations
        KGeoBag::KThreeVector E=fCharge / (KConst::FourPiEps() * c *pow(1.- fReceiverDir.Dot(fNewVelocity)/c,3.)*fReceiverDistance ) * (c / (fGamma*fGamma*fReceiverDistance)*(fReceiverDir-1./c*fNewVelocity) + (fReceiverDir.Cross(fReceiverDir.Cross(BetaDot))));

        KGeoBag::KThreeVector B = fReceiverDir.Cross(E)/c;

        KGeoBag::KThreeVector S =1. / KConst::MuNull()*E.Cross(B);

        //double P=fabs(S.Dot(SurfaceNormal));
        double P=fabs(S.Dot(SurfaceNormal))/fabs(fCharge);

        return P;
    }

    void ParticleSlim::Print()
    {
        printf("-----------------------------\n");
        printf("fTime: %e\n",fTime);
        printf("fTimeStep: %e\n",fTimeStep);

        printf("fPosition: %e %e %e\n",fPosition[0],fPosition[1],fPosition[2]);
        printf("fVelocity: %e %e %e\n",fVelocity[0],fVelocity[1],fVelocity[2]);
        printf("fVelocityParallel: %e\n",fVelocityParallel);

        printf("fNewPosition: %e %e %e\n",fNewPosition[0],fNewPosition[1],fNewPosition[2]);
        printf("fNewVelocity: %e %e %e\n",fNewVelocity[0],fNewVelocity[1],fNewVelocity[2]);

        printf("fGuidingCenterPosition: %e %e %e\n",fGuidingCenterPosition[0],fGuidingCenterPosition[1],fGuidingCenterPosition[2]);
        printf("fMagneticField: %e %e %e\n",fMagneticField[0],fMagneticField[1],fMagneticField[2]);

        printf("fAlpha: %e %e %e\n",fAlpha[0],fAlpha[1],fAlpha[2]);
        printf("fBeta: %e %e %e\n",fBeta[0],fBeta[1],fBeta[2]);

        printf("fMass: %e\n",fMass);
        printf("fCharge: %e\n",fCharge);
        printf("fGamma: %e\n",fGamma);
        printf("fCyclotronFrequency: %e\n",fCyclotronFrequency);
        printf("fCyclotronRadius: %e\n",fCyclotronRadius);

        /////////////////Receiver Info////////////
        printf("fReceiverTime: %e\n",fReceiverTime);
        printf("fReceiverPosition: %e %e %e\n",fReceiverPosition[0],fReceiverPosition[1],fReceiverPosition[2]);
        printf("fReceiverVector: %e %e %e\n",fReceiverVector[0],fReceiverVector[1],fReceiverVector[2]);
        printf("fReceiverDistance: %e\n",fReceiverDistance);
        printf("-----------------------------\n");
        
        return;
    }

    void ParticleSlim::SetSpline(ParticleSlim aNextParticle)
    {
        double dt=aNextParticle.GetTime()-fTime;
        fTimeStep=dt;

        double CyclotronFrequency=fCyclotronFrequency;
        if(fCharge<0)CyclotronFrequency*=-1.;

        fSplineD= dt * (aNextParticle.GetVelocity()-fVelocityParallel*fMagneticField.Unit()-fCyclotronRadius*fCyclotronFrequency*(-sin(CyclotronFrequency*dt)*fAlpha+cos(CyclotronFrequency*dt)*fBeta)) - 3.*(aNextParticle.GetPosition()-fGuidingCenterPosition-fVelocityParallel*fMagneticField.Unit()*dt-fCyclotronRadius*(cos(CyclotronFrequency*dt)*fAlpha+sin(CyclotronFrequency*dt)*fBeta));
        fSplineC= aNextParticle.GetPosition()-fGuidingCenterPosition-fVelocityParallel*fMagneticField.Unit()*dt-fCyclotronRadius*(cos(CyclotronFrequency*dt)*fAlpha+sin(CyclotronFrequency*dt)*fBeta)-fSplineD;

        //KGeoBag::KThreeVector newPos = fGuidingCenterPosition + fVelocityParallel * fMagneticField.Unit() * dt + fCyclotronRadius * ( cos(CyclotronFrequency*dt)*fAlpha + sin(CyclotronFrequency*dt)*fBeta) + fSplineC + fSplineD;

        //KGeoBag::KThreeVector newVel=fVelocityParallel*fMagneticField.Unit() + fCyclotronFrequency*fCyclotronRadius*(-sin(CyclotronFrequency*dt)*fAlpha+cos(CyclotronFrequency*dt)*fBeta)+3.*fSplineC/dt+4.*fSplineD/dt;
        //printf("%e %e\n",fSplineC.Magnitude(),fSplineD.Magnitude());

        //KGeoBag::KThreeVector tmp=newPos-aNextParticle.GetPosition();
        //KGeoBag::KThreeVector tmp=newVel-aNextParticle.GetVelocity();
        //printf("%e\n",tmp.Magnitude());

        //fSplineC/=(dt*dt*dt*dt);
        //fSplineD/=(dt*dt*dt);

        return;
    }

    void ParticleSlim::ErrorCheck(ParticleSlim a)
    {
        //KGeoBag::KThreeVector tmp=fNewVelocity-a.fVelocity;
        KGeoBag::KThreeVector tmp=fNewPosition-a.fPosition;
        printf("%e\n",tmp.Magnitude());

    }

} /* namespace locust */
