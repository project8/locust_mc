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
            void SetPosition(double,double,double);
            void SetVelocityVector(double,double,double);
            void SetMagneticFieldVector(double,double,double);
            void SetMass(double);
            void SetCharge(double);
            void SetCyclotronFrequency(double);

            void SetReceiverPosition(double,double,double);
            void SetReceiverTime(double);
            void SetKinematicProperties();

            void Interpolate(double);

            double GetTime();
            double GetTimeDisplacement();
            KGeoBag::KThreeVector GetPosition();
            KGeoBag::KThreeVector GetVelocity();

            double GetSpaceTimeInterval();
            double GetSpaceTimeInterval(double);
            double GetReceiverDistance();

            double CalculatePower(double,double,double);
            double CalculateVoltage();

            void Print();
            void SetSpline(ParticleSlim);
            void ErrorCheck(ParticleSlim);



        //private:
            double fTime; //Time of node given from kassiopeia
            double fTimeDisplacement;//Allows us to interpolate particle locally around node
            double fTimeStep;//Time Step to next node

            KGeoBag::KThreeVector fPosition;
            KGeoBag::KThreeVector fVelocity;
            double fVelocityParallel;

            KGeoBag::KThreeVector fNewPosition;
            KGeoBag::KThreeVector fNewVelocity;

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

            KGeoBag::KThreeVector fSplineC;
            KGeoBag::KThreeVector fSplineD;

            /////////////////Receiver Info////////////
            double fReceiverTime;
            KGeoBag::KThreeVector fReceiverPosition;
            KGeoBag::KThreeVector fReceiverVector;
            double fReceiverDistance;
            KGeoBag::KThreeVector fReceiverDir;

};


} /* namespace locust */

#endif /* LMCPARTICLESLIM_HH_ */
