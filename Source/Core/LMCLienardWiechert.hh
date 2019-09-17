/*
 * LMCLienardWiechert.hh
 *
 *  Created on: Sep 17, 2019
 *      Author: nbuzinsky
 */

#ifndef LMCLIENARDWIECHERT_HH_
#define LMCLIENARDWIECHERT_HH_

#include "LMCParticle.hh"
#include "LMCThreeVector.hh"
#include <vector>

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
            void SetFieldEvent(const double aTime, const LMCThreeVector aFieldPoint);
            LMCThreeVector GetElectricField();
            LMCThreeVector GetMagneticField();


        private:
            int FindNode(double tNew) const;
            double GetSpaceTimeInterval(const double &aParticleTime, const double &aReceiverTime, const LMCThreeVector &aParticlePosition, const LMCThreeVector &aReceiverPosition ) const;
            double GetStepRoot(const locust::Particle aParticle, double aReceiverTime, LMCThreeVector aReceiverPosition, double aSpaceTimeInterval) const;
            void CacheSolution(const int aCurrentIndex, const double aRetardedTime);
            void InitialRetardedTimeGuess();
            bool IsInLightCone();

            std::vector<std::pair<int, double> > fPreviousTimes; //Cache the results from previous iteration. [0] is previous index, [1] is corresponding retarded time of previous solution
            LMCParticle fCurrentParticle;
            LMCThreeVector fFieldPosition;
            double fFieldTime;
            unsigned fPatchIndex;

            ///////////////

};


} /* namespace locust */

#endif /* LMCLIENARDWIECHERT_HH_ */
