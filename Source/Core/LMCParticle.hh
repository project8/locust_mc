/*
 * LMCParticle.hh
 *
 *  Created on: Jan 13, 2017
 *      Author: nbuzinsky
 */

#ifndef LMCPARTICLE_HH_
#define LMCPARTICLE_HH_

#include <stdint.h>
#include <KThreeVector.hh>
#include <KConst.h>


namespace locust
{

    class Particle
    {

        public:
            Particle();
            virtual ~Particle();

            void SetTime(double);
            void SetPosition(double,double,double);
            void SetVelocityVector(double,double,double);
            void SetMagneticFieldVector(double,double,double);
            void SetCyclotronFrequency(double);
            void SetMass(double);
            void SetCharge(double);

            void SetKinematicProperties();

            void Interpolate(double);

            double GetTime(const bool &aInterpolated = false) const;
            KGeoBag::KThreeVector GetPosition(const bool &aInterpolated = false ) const;
            KGeoBag::KThreeVector GetVelocity(const bool &aInterpolated = false) const;
            KGeoBag::KThreeVector GetAcceleration(const bool &aInterpolated = false) const;
            double GetCyclotronFrequency() const;
            double GetKineticEnergy() const;
            double GetCharge() const;
            double GetLarmorPower() const;

            double CalculateVoltage(const KGeoBag::KThreeVector&) const;
            KGeoBag::KThreeVector CalculateElectricField(const KGeoBag::KThreeVector&) const;
            KGeoBag::KThreeVector CalculateMagneticField(const KGeoBag::KThreeVector&) const;

            void Print() const;
            void SetSpline(Particle);

        private:
            double fTime; //Time of node given from kassiopeia
            double fTimeDisplacement;//Allows us to interpolate particle locally around node
            double fTimeStep;//Time Step to next node

            KGeoBag::KThreeVector fPosition;
            KGeoBag::KThreeVector fVelocity;
            KGeoBag::KThreeVector fAcceleration;
            double fVelocityParallel;

            KGeoBag::KThreeVector fNewPosition;
            KGeoBag::KThreeVector fNewVelocity;
            KGeoBag::KThreeVector fNewAcceleration;

            KGeoBag::KThreeVector fGuidingCenterPosition;
            KGeoBag::KThreeVector fMagneticField;

            //2 perp. vectors which define helical motion
            KGeoBag::KThreeVector fAlpha;
            KGeoBag::KThreeVector fBeta;

            double fMass;
            double fCharge;
            double fGamma;
            double fCyclotronFrequency;
            double fCyclotronRadius;
            double fLarmorPower;

            KGeoBag::KThreeVector fSplineC;
            KGeoBag::KThreeVector fSplineD;

};


} /* namespace locust */

#endif /* LMCPARTICLE_HH_ */
