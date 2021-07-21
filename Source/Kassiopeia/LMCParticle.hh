/*
 * LMCParticle.hh
 *
 *  Created on: Jan 13, 2017
 *      Author: nbuzinsky
 */

#ifndef LMCPARTICLE_HH_
#define LMCPARTICLE_HH_

#include "LMCConst.hh"
#include <stdint.h>
#include "LMCThreeVector.hh"


namespace locust
{
 /*!
 @class Particle
 @author N. Buzinsky
 @brief Used to describe the properties of a particle
 @details
 Available configuration options:
 No input parameters
 */
    class Particle
    {

        public:
            Particle();
            virtual ~Particle();

            void SetTime(double);
            void SetPosition(double,double,double);
            void SetVelocityVector(double,double,double);
            void SetMagneticFieldVector(double,double,double);
            void SetPitchAngle(double);
            void SetCyclotronFrequency(double);
            void SetMass(double);
            void SetCharge(double);

            void SetKinematicProperties();

            void Interpolate(double);

            double GetTime(const bool &aInterpolated = false) const;
            LMCThreeVector GetPosition(const bool &aInterpolated = false ) const;
            LMCThreeVector GetVelocity(const bool &aInterpolated = false) const;
            LMCThreeVector GetAcceleration(const bool &aInterpolated = false) const;

            double GetPitchAngle() const;
            double GetCyclotronFrequency() const;
            double GetKineticEnergy() const;
            double GetCharge() const;
            double GetLarmorPower() const;

            double CalculateVoltage(const LMCThreeVector&) const;
            LMCThreeVector CalculateElectricField(const LMCThreeVector&) const;
            LMCThreeVector CalculateMagneticField(const LMCThreeVector&) const;

            void Print() const;
            void SetSpline(Particle);

        private:
            double fTime; //Time of node given from kassiopeia
            double fTimeDisplacement;//Allows us to interpolate particle locally around node
            double fTimeStep;//Time Step to next node

            LMCThreeVector fPosition;
            LMCThreeVector fVelocity;
            LMCThreeVector fAcceleration;
            double fVelocityParallel;

            LMCThreeVector fNewPosition;
            LMCThreeVector fNewVelocity;
            LMCThreeVector fNewAcceleration;

            LMCThreeVector fGuidingCenterPosition;
            LMCThreeVector fMagneticField;

            double fPitchAngle;

            //2 perp. vectors which define helical motion
            LMCThreeVector fAlpha;
            LMCThreeVector fBeta;

            double fMass;
            double fCharge;
            double fGamma;
            double fCyclotronFrequency;
            double fCyclotronRadius;
            double fLarmorPower;

            LMCThreeVector fSplineC;
            LMCThreeVector fSplineD;

};


} /* namespace locust */

#endif /* LMCPARTICLE_HH_ */
