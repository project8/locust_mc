/*
 * LMCParticleSlim.hh
 *
 *  Created on: Jan 13, 2017
 *      Author: nbuzinsky
 */

#ifndef LMCPARTICLESLIM_HH_
#define LMCPARTICLESLIM_HH_

#include <stdint.h>
#include <KThreeVector.hh>
#include <KConst.h>


namespace locust
{

    class ParticleSlim
    {

        public:
            ParticleSlim();
            virtual ~ParticleSlim();

            void SetTime(double);
            double GetTime();
            void SetPosition(double,double,double);
            void SetVelocityVector(double,double,double);
            void SetMagneticFieldVector(double,double,double);
            void SetMass(double);
            void SetCharge(double);

            void SetReceiverTime(double);
            void SetReceiverPosition(double,double,double);
            void CalculateKinematicProperties();

            void CalculateQuadraticCoefficients(double&,double&,double&);
            double Interpolate(double);
            double GetSpaceTimeInterval();
            double GetSpaceTimeInterval(double,KGeoBag::KThreeVector);
            double GetReceiverDistance();
            double NewtonStep(double, double);
            double HouseHolderStep(double, double);

            double CalculatePower(double,double,double);




        private:
            double fTime;

            KGeoBag::KThreeVector fPosition;
            KGeoBag::KThreeVector fVelocity;
            KGeoBag::KThreeVector fBetaVelocity;
            double fVelocityParallel;

            KGeoBag::KThreeVector fGuidingCenterPosition;

            KGeoBag::KThreeVector fMagneticField;

            KGeoBag::KThreeVector fAlpha;
            KGeoBag::KThreeVector fBeta;

            double fMass;
            double fCharge;
            double fGamma;
            double fCyclotronFrequency;
            double fCyclotronRadius;

            /////////////////Receiver Info////////////
            double fReceiverTime;
            KGeoBag::KThreeVector fReceiverPosition;
            KGeoBag::KThreeVector fReceiverVector;
            double fReceiverDistance;
            KGeoBag::KThreeVector fReceiverDir;

};


} /* namespace locust */

#endif /* LMCPARTICLESLIM_HH_ */
