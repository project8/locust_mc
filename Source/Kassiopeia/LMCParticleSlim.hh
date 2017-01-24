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

            void SetReceiverPosition(double,double,double);
            void SetReceiverTime(double);
            void SetKinematicProperties();

            void Interpolate(double);

            void CalculateQuadraticCoefficients(double&,double&,double&);
            double GetSpaceTimeInterval();
            double GetSpaceTimeInterval(double);
            double GetReceiverDistance();
            double NewtonStep(double);
            double HouseHolderStep(double);

            double CalculatePower(double,double,double);

            void Print();



        private:
            double fTime; //Time od node given from kassiopeia
            double fTimeDisplacement;//dt. Allows us to interpolate particle locally around node

            KGeoBag::KThreeVector fPosition;
            KGeoBag::KThreeVector fVelocity;
            double fVelocityParallel;

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

            /////////////////Receiver Info////////////
            double fReceiverTime;
            KGeoBag::KThreeVector fReceiverPosition;
            KGeoBag::KThreeVector fReceiverVector;
            double fReceiverDistance;
            KGeoBag::KThreeVector fReceiverDir;

};


} /* namespace locust */

#endif /* LMCPARTICLESLIM_HH_ */
