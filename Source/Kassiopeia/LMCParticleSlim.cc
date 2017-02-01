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
            fPosition( -99., -99., -99. ),
            fVelocity( -99., -99., -99. ),
            fVelocityParallel( -99.),
            fGuidingCenterPosition( -99., -99., -99. ),
            fMagneticField( -99., -99., -99. ),
            fAlpha( -99., -99., -99. ),
            fBeta( -99., -99., -99. ),
            fMass( -99. ),
            fCharge( -99. ),
            fGamma( -99. ),
            fCyclotronFrequency( -99. ),
            fCyclotronRadius( -99. ),
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
        return;
    }

    void ParticleSlim::SetVelocityVector(double Vx, double Vy, double Vz)
    {
        fVelocity.SetComponents(Vx,Vy,Vz);
        return;
    }

    void ParticleSlim::SetMagneticFieldVector(double Bx, double By, double Bz)
    {
        fMagneticField.SetComponents(Bx,By,Bz);
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
        fReceiverVector=fReceiverPosition-fPosition;
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
        fCyclotronFrequency=fabs(fCharge) * fMagneticField.Magnitude() / ( fMass * fGamma);

        fVelocityParallel=fVelocity.Dot(fMagneticField.Unit());

        KGeoBag::KThreeVector vPerp=fVelocity-fVelocityParallel*fMagneticField.Unit();
        fCyclotronRadius = fGamma * fMass * vPerp.Magnitude() / (fabs(fCharge) * fMagneticField.Magnitude());

        KGeoBag::KThreeVector tCrossRadial=(fVelocity.Cross(fMagneticField)).Unit();

        if(fCharge<0)
        {
            tCrossRadial=-1.*tCrossRadial;
        }

        fGuidingCenterPosition=fPosition+fCyclotronRadius*tCrossRadial;

        fAlpha=fPosition-fGuidingCenterPosition;

        fBeta=fMagneticField.Cross(fAlpha);
        fAlpha=fAlpha.Unit();
        fBeta=fBeta.Unit();
        return;
    }

    void ParticleSlim::Interpolate(double tNew)
    {
        double dt=tNew-fTime;
        fTimeDisplacement=dt;
        fPosition=fGuidingCenterPosition + fCyclotronRadius * ( cos(fCyclotronFrequency*fTimeDisplacement)*fAlpha + sin(fCyclotronFrequency*fTimeDisplacement)*fBeta);
        fVelocity=fVelocityParallel*fMagneticField.Unit()+ fCyclotronRadius*fCyclotronFrequency*(-sin(fCyclotronFrequency*fTimeDisplacement)*fAlpha+cos(fCyclotronFrequency*fTimeDisplacement)*fBeta);

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
        KGeoBag::KThreeVector newPos = fGuidingCenterPosition + fCyclotronRadius * ( cos(fCyclotronFrequency*dt)*fAlpha + sin(fCyclotronFrequency*dt)*fBeta);
        KGeoBag::KThreeVector dr=fReceiverPosition-newPos;
        double s=fReceiverTime-aTime-dr.Magnitude()/KConst::C();
        return s;
    }


    double ParticleSlim::GetReceiverDistance()
    {
        return fReceiverDistance;
    }


    double ParticleSlim::CalculatePower(double X, double Y, double Z)
    {
        double c=KConst::C();

        KGeoBag::KThreeVector SurfaceNormal(X,Y,Z);

        //Lorentz Equation
        KGeoBag::KThreeVector BetaDot=fCharge/(fMass * c) *fVelocity.Cross(fMagneticField);

        //Lienard-Wiechert Equations
        KGeoBag::KThreeVector E=fCharge / (KConst::FourPiEps() * c *pow(1.- fReceiverDir.Dot(fVelocity)/c,3.)*fReceiverDistance ) * (c / (fGamma*fGamma*fReceiverDistance)*(fReceiverDir-1./c*fVelocity) + (fReceiverDir.Cross(fReceiverDir.Cross(BetaDot))));

        KGeoBag::KThreeVector B = fReceiverDir.Cross(E)/c;

        KGeoBag::KThreeVector S =1. / KConst::MuNull()*E.Cross(B);

        double P=S.Dot(SurfaceNormal);

        return P;
    }

    void ParticleSlim::Print()
    {
        printf("-----------------------------\n");
        printf("fTime: %e\n",fTime);

        printf("fPosition: %e %e %e\n",fPosition[0],fPosition[1],fPosition[2]);
        printf("fVelocity: %e %e %e\n",fVelocity[0],fVelocity[1],fVelocity[2]);
        printf("fVelocityParallel: %e\n",fVelocityParallel);

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


} /* namespace locust */
