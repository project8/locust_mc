/*
 * LMCParticle.cc
 *
 *  Created on: Jan 13, 2017
 *      Author: nbuzinsky
 */

#include "LMCParticle.hh"


namespace locust
{

    Particle::Particle() :
            fTime( -99. ), 
            fTimeDisplacement( -99.), 
            fTimeStep( -99.), 
            fPosition( -99., -99., -99. ),
            fVelocity( -99., -99., -99. ),
            fAcceleration( -99., -99., -99. ),
            fVelocityParallel( -99.),
            fNewPosition( -99., -99., -99. ),
            fNewVelocity( -99., -99., -99. ),
            fNewAcceleration( -99., -99., -99. ),
            fGuidingCenterPosition( -99., -99., -99. ),
            fMagneticField( -99., -99., -99. ),
            fPitchAngle( -99. ),
            fAlpha( -99., -99., -99. ),
            fBeta( -99., -99., -99. ),
            fMass( -99. ),
            fCharge( -99. ),
            fGamma( -99. ),
            fCyclotronFrequency( -99. ),
            fCyclotronRadius( -99. ),
            fLarmorPower( -99. ),
            fSplineC( -99., -99., -99. ),
            fSplineD( -99., -99., -99. )
    {
    }

    Particle::~Particle()
    {

    }

    void Particle::SetTime(double aTime)
    {
        fTime = aTime;
        fTimeDisplacement = 0.;
        return;
    }
    double Particle::GetTime(const bool &aInterpolated) const
    {
        return aInterpolated ? fTime + fTimeDisplacement : fTime;
    }

    void Particle::SetPosition(double X, double Y, double Z)
    {
        fPosition.SetComponents(X,Y,Z);
        fNewPosition.SetComponents(X,Y,Z);
        return;
    }

    LMCThreeVector Particle::GetPosition(const bool &aInterpolated ) const
    {
        return aInterpolated ? fNewPosition : fPosition;
    }

    void Particle::SetVelocityVector(double Vx, double Vy, double Vz)
    {
        fVelocity.SetComponents(Vx,Vy,Vz);
        fNewVelocity.SetComponents(Vx,Vy,Vz);
        return;
    }

    LMCThreeVector Particle::GetVelocity(const bool &aInterpolated ) const
    {
        return aInterpolated ? fNewVelocity : fVelocity;
    }

    LMCThreeVector Particle::GetAcceleration(const bool &aInterpolated ) const
    {
        return aInterpolated ? fNewAcceleration : fNewAcceleration;
    }


    void Particle::SetMagneticFieldVector(double Bx, double By, double Bz)
    {
        fMagneticField.SetComponents(Bx,By,Bz);
        return;
    }
    void Particle::SetCyclotronFrequency(double aCyclotronFrequency)
    {
        fCyclotronFrequency=aCyclotronFrequency;
        return;
    }

    double Particle::GetCyclotronFrequency() const
    {
        return fCyclotronFrequency;
    }

    void Particle::SetPitchAngle(double aPitchAngle)
    {
        fPitchAngle=aPitchAngle;
        return;
    }



    double Particle::GetPitchAngle() const
    {
        return fPitchAngle;
    }


    void Particle::SetMass(double aMass)
    {
        fMass = aMass;
        return;
    }

    void Particle::SetCharge(double aCharge)
    {
        fCharge = aCharge;
        return;
    }

    double Particle::GetCharge() const
    {
        return fCharge;
    }

    double Particle::GetLarmorPower() const
    {
        return fLarmorPower;
    }

    void Particle::SetKinematicProperties()
    {
        fAcceleration = fCharge/fMass*fNewVelocity.Cross(fMagneticField);
        fNewAcceleration = fAcceleration;

        fLarmorPower = 2. / 3. * pow(fCharge , 2.) * fAcceleration.MagnitudeSquared() / (LMCConst::FourPiEps() * pow( LMCConst::C() , 3.) );

        double tBeta = fVelocity.Magnitude() / LMCConst::C();
        fGamma = 1. / sqrt(( 1. - tBeta ) * ( 1. + tBeta ));

        fVelocityParallel=fVelocity.Dot(fMagneticField.Unit());
        LMCThreeVector vPerp = fVelocity - fVelocityParallel * fMagneticField.Unit();
        fBeta=vPerp.Unit();

        fCyclotronRadius = fGamma * fMass * vPerp.Magnitude() / (fabs(fCharge) * fMagneticField.Magnitude());

        LMCThreeVector tCrossRadial = (fVelocity.Cross(fMagneticField)).Unit();

        if(fCharge<0) tCrossRadial *= -1.;

        fGuidingCenterPosition = fPosition + fCyclotronRadius * tCrossRadial;

        fAlpha = fPosition - fGuidingCenterPosition;

        fAlpha=fAlpha.Unit();
        //fBeta=vPerp.Unit();

        fSplineC.SetComponents(0.,0.,0.);
        fSplineD.SetComponents(0.,0.,0.);

        return;
    }

    void Particle::Interpolate(double tNew)
    {
        fTimeDisplacement = tNew - fTime; 
        
        double tCyclotronFrequency = fCyclotronFrequency;
        if(fCharge<0) tCyclotronFrequency *= -1.;

        fNewPosition = fGuidingCenterPosition + fVelocityParallel * fMagneticField.Unit() * fTimeDisplacement + fCyclotronRadius * ( cos( tCyclotronFrequency * fTimeDisplacement) * fAlpha + sin( tCyclotronFrequency * fTimeDisplacement) * fBeta) + fSplineC * pow( fTimeDisplacement / fTimeStep , 3. ) + fSplineD * pow( fTimeDisplacement / fTimeStep, 4. );

        fNewVelocity = fVelocityParallel * fMagneticField.Unit() + fCyclotronRadius * fCyclotronFrequency *( - sin(tCyclotronFrequency * fTimeDisplacement) * fAlpha + cos( tCyclotronFrequency * fTimeDisplacement ) * fBeta ) + 3. * fSplineC * pow( fTimeDisplacement / fTimeStep, 2. ) / fTimeStep + 4. * fSplineD * pow( fTimeDisplacement / fTimeStep, 3. ) / fTimeStep;

        fNewAcceleration = fCharge / fMass * fNewVelocity.Cross(fMagneticField);

        return;
    }


    LMCThreeVector Particle::CalculateElectricField(const LMCThreeVector &aFieldPosition) const
    {
        double c=LMCConst::C();

        LMCThreeVector tFieldPositionVector =  aFieldPosition - fNewPosition;
        LMCThreeVector tFieldPositionNormal =  tFieldPositionVector.Unit();
        double tFieldPositionDistance = tFieldPositionVector.Magnitude();
        //Lorentz Equation
        LMCThreeVector tBetaDot = fNewAcceleration / c;
        LMCThreeVector tNormalMinusBeta = tFieldPositionNormal - 1./ c * fNewVelocity;

        //Lienard-Wiechert Equations
        LMCThreeVector E = fCharge * (tFieldPositionNormal.Cross(tNormalMinusBeta.Cross(tBetaDot) ) ) / ( LMCConst::FourPiEps() * c *pow( 1. - tFieldPositionNormal.Dot(fNewVelocity) / c, 3. ) * tFieldPositionDistance );
        return E;
    }

    LMCThreeVector Particle::CalculateMagneticField(const LMCThreeVector &aFieldPosition) const
    {
        LMCThreeVector tFieldPositionVector =  aFieldPosition - fNewPosition;
        LMCThreeVector tFieldPositionNormal =  tFieldPositionVector.Unit();
        
        LMCThreeVector H = tFieldPositionNormal.Cross(CalculateElectricField(aFieldPosition)) / ( LMCConst::MuNull() * LMCConst::C() );

        return H;
    }

    double Particle::CalculateVoltage(const LMCThreeVector &aFieldPosition) const
    {
        //Lienard-Wiechert Equations
        LMCThreeVector tFieldPositionVector =  aFieldPosition - fNewPosition;
        LMCThreeVector tFieldPositionNormal =  tFieldPositionVector.Unit();
        double tFieldPositionDistance = tFieldPositionVector.Magnitude();

        double V=fCharge / (LMCConst::FourPiEps()*tFieldPositionDistance*(1.- tFieldPositionNormal.Dot(fNewVelocity)/LMCConst::C()));
        return V;

    }

    void Particle::Print() const
    {
        printf("-----------------------------------------------------\n");
        printf("fTime: %e\n",fTime);
        printf("fTimeStep: %e\n",fTimeStep);

        printf("fPosition: %e %e %e\n",fPosition[0],fPosition[1],fPosition[2]);
        printf("fVelocity: %e %e %e\n",fVelocity[0],fVelocity[1],fVelocity[2]);
        printf("fVelocityParallel: %e\n",fVelocityParallel);
        printf("fAcceleration: %e %e %e\n",fAcceleration[0],fAcceleration[1],fAcceleration[2]);

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
        printf("-----------------------------------------------------\n");

        
        return;
    }

    void Particle::SetSpline(Particle aNextParticle)
    {
        double dt=aNextParticle.GetTime()-fTime;
        fTimeStep=dt;

        double CyclotronFrequency=fCyclotronFrequency;
        if(fCharge<0)CyclotronFrequency*=-1.;

        fSplineD= dt * (aNextParticle.GetVelocity()-fVelocityParallel*fMagneticField.Unit()-fCyclotronRadius*fCyclotronFrequency*(-sin(CyclotronFrequency*dt)*fAlpha+cos(CyclotronFrequency*dt)*fBeta)) - 3.*(aNextParticle.GetPosition()-fGuidingCenterPosition-fVelocityParallel*fMagneticField.Unit()*dt-fCyclotronRadius*(cos(CyclotronFrequency*dt)*fAlpha+sin(CyclotronFrequency*dt)*fBeta));
        fSplineC= aNextParticle.GetPosition()-fGuidingCenterPosition-fVelocityParallel*fMagneticField.Unit()*dt-fCyclotronRadius*(cos(CyclotronFrequency*dt)*fAlpha+sin(CyclotronFrequency*dt)*fBeta)-fSplineD;

        return;
    }

} /* namespace locust */
