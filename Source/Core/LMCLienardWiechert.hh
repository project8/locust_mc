/*
 * LMCLienardWiechert.hh
 *
 *  Created on: Sep 17, 2019
 *      Author: nbuzinsky
 */

#ifndef LMCLIENARDWIECHERT_HH_
#define LMCLIENARDWIECHERT_HH_

#include "LMCThreeVector.hh"

namespace locust
{
 /*!
 @class LienardWiechert
 @author N. Buzinsky
 @brief Helper class that computes fields at arbitrary point (t,x,y,z) given particle trajectory (from fParticleHistory)
 @details
 Available configuration options:
 No input parameters
 */
    class LienardWeichert
    {

        public:
            LienardWiechert();
            virtual ~LienardWiechert();
            void SolveFieldSolutions();
            void SetFieldPoint(double aTime, LMCThreeVector aFieldPoint);
            LMCThreeVector GetElectricField();
            LMCThreeVector GetMagneticField();
            //void SetFieldPoint(double aTime, double aFieldPointX, double aFieldPointY, double aFieldPointZ);


        private:
            int FindNode(double tNew) const;
            double GetSpaceTimeInterval(const double &aParticleTime, const double &aReceiverTime, const LMCThreeVector &aParticlePosition, const LMCThreeVector &aReceiverPosition ) const;
            double GetStepRoot(const locust::Particle aParticle, double aReceiverTime, LMCThreeVector aReceiverPosition, double aSpaceTimeInterval) const;
            void CacheSolution();
            void InitialRetardedTimeGuess();
            bool IsInLightCone();

        const double kassiopeiaTimeStep = fabs(fParticleHistory[0].GetTime() - fParticleHistory[1].GetTime());
        const int historySize = fParticleHistory.size();
        locust::Particle tCurrentParticle = fParticleHistory.back();
        int CurrentIndex;
        double tSpaceTimeInterval=99.;
        double dtRetarded=0;
        double tTolerance=1e-23;


};


} /* namespace locust */

#endif /* LMCPARTICLE_HH_ */
