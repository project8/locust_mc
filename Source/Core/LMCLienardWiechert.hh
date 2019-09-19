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
    class LienardWiechert
    {

        public:
            LienardWiechert();
            virtual ~LienardWiechert();

            void AddFieldPoint(const LMCThreeVector aFieldPoint);
            void SetFieldEvent(const double aTime, const unsigned aFieldPointIndex);
            boolt SolveFieldSolutions();

            LMCThreeVector GetElectricField();
            LMCThreeVector GetMagneticField();

        private:
            bool IsInLightCone() const;
            std::pair<unsigned, double> GuessRetardedTime()  const;
            double GetSpaceTimeInterval(const double &aParticleTime, const double &aReceiverTime, const LMCThreeVector &aParticlePosition, const LMCThreeVector &aReceiverPosition ) const;
            double GetStepRoot(const locust::Particle aParticle, double aReceiverTime, LMCThreeVector aReceiverPosition, double aSpaceTimeInterval) const;
            unsigned FindClosestParticle(double tNew) const;
            std::pair<unsigned, double> FindRoot() const;
            void CacheSolution(const int aCurrentIndex, const double aRetardedTime);

            LMCParticle fCurrentParticle;
            LMCThreeVector fFieldPosition;
            double fFieldTime;
            unsigned fAntennaIndex;

            std::vector<LMCThreeVector> fAntennaPositions;
            std::vector<std::pair<unsigned, double> > fCachedSolutions; //Cache the results from previous iteration. [0] is previous index, [1] is corresponding retarded time of previous solution
            std::vector<bool> fHasCachedSolution;

};


} /* namespace locust */

#endif /* LMCLIENARDWIECHERT_HH_ */
