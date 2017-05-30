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
            fNewAcceleration( -99., -99., -99. ),
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
        fNewAcceleration=fCharge/fMass*fNewVelocity.Cross(fMagneticField);

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

        fNewAcceleration=fCharge/fMass*fNewVelocity.Cross(fMagneticField);

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

    double ParticleSlim::GetStepRoot(int StepOrder, double tSpaceTimeInterval)
    {
        double tRetarded=fTime+fTimeDisplacement;

        double c=KConst::C();


        if(StepOrder==0)
        {
            return tRetarded+tSpaceTimeInterval;
        }

        double fZero=pow((fReceiverTime-tRetarded),2.)-pow(fReceiverDistance,2.)/(c*c);
        double fZeroPrime=2.*((tRetarded-fReceiverTime)-fNewVelocity.Dot(fNewPosition)/(c*c)+fNewVelocity.Dot(fReceiverPosition)/(c*c));
        double NewtonRatio=fZero/fZeroPrime;

        if(StepOrder==1)
        {
            return tRetarded-NewtonRatio;
        }

        double fZeroDoublePrime=2.*(1.-fNewVelocity.Dot(fNewVelocity)/(c*c)-fNewAcceleration.Dot(fNewPosition-fReceiverPosition)/(c*c));

        if(StepOrder==2)
        {
            return tRetarded-NewtonRatio*(1.+(NewtonRatio*fZeroDoublePrime)/(2.*fZeroPrime));
        }

        return 0;
    }


    void ParticleSlim::CalculateElectricField(double& Ex, double &Ey, double& Ez)
    {
        double c=KConst::C();

        //Lorentz Equation
        KGeoBag::KThreeVector BetaDot = fNewAcceleration / c;
        KGeoBag::KThreeVector nMinusBeta = fReceiverDir - 1./ c * fNewVelocity;

        //Lienard-Wiechert Equations
        KGeoBag::KThreeVector E=fCharge / (KConst::FourPiEps() * c *pow(1.- fReceiverDir.Dot(fNewVelocity)/c,3.)*fReceiverDistance ) * ( c * nMinusBeta / (fGamma*fGamma*fReceiverDistance) + (fReceiverDir.Cross(nMinusBeta.Cross(BetaDot))));

        Ex=E.X();
        Ey=E.Y();
        Ez=E.Z();


    }

    void ParticleSlim::CalculateMagneticField(double& Hx, double &Hy, double& Hz)
    {
        double Ex, Ey, Ez;
        CalculateElectricField(Ex,Ey,Ez);
        KGeoBag::KThreeVector E(Ex,Ey,Ez);
        
        KGeoBag::KThreeVector H = fReceiverDir.Cross(E)/(KConst::MuNull()*KConst::C());

        Hx=H.X();
        Hy=H.Y();
        Hz=H.Z();
    }

    double ParticleSlim::CalculateVoltage()
    {
        //Lienard-Wiechert Equations
        double V=fCharge / (KConst::FourPiEps()*fReceiverDistance*(1.- fReceiverDir.Dot(fNewVelocity)/KConst::C()));
        //double V=-1./ (fReceiverDistance*(1.- fReceiverDir.Dot(fNewVelocity)/KConst::C()));
        return V;

    }

    ////Only works with sphere at origin
    double ParticleSlim::CalculatePower()
    {
        double Ex, Ey, Ez;
        CalculateElectricField(Ex,Ey,Ez);
        KGeoBag::KThreeVector E(Ex,Ey,Ez);
        KGeoBag::KThreeVector H = fReceiverDir.Cross(E)/(KConst::MuNull()*KConst::C());
        KGeoBag::KThreeVector S = E.Cross(H);

        return S.Dot(fReceiverPosition.Unit())*(1.-fReceiverDir.Dot(fNewVelocity)/KConst::C());
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
        printf("fNewAcceleration: %e %e %e\n",fNewAcceleration[0],fNewAcceleration[1],fNewAcceleration[2]);

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

        return;
    }

    void ParticleSlim::ErrorCheck(ParticleSlim a)
    {
        //KGeoBag::KThreeVector tmp=fNewVelocity-a.fVelocity;
        KGeoBag::KThreeVector tmp=fNewPosition-a.fPosition;
        printf("%e\n",tmp.Magnitude());

    }

} /* namespace locust */
