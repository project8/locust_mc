/*
 * LMCLienardWiechert.hh
 *
 *  Created on: Sep 17, 2019
 *      Author: nbuzinsky
 */

#ifndef LMCLIENARDWIECHERT_HH_
#define LMCLIENARDWIECHERT_HH_

#include "LMCKassLocustInterface.hh"
#include "LMCParticle.hh"
#include "LMCThreeVector.hh"
#include <algorithm>
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

            void SetParticleHistory( std::deque<Particle> aParticleHistory );
            void AddFieldPoint(const LMCThreeVector aFieldPoint);
            void SetFieldEvent(const double aTime, const unsigned aFieldPointIndex);
            void SolveFieldSolutions();

            LMCThreeVector GetElectricField( ) const;
            LMCThreeVector GetMagneticField() const;
            locust::Particle GetRetardedParticle() const;

        private:
            void SetKassiopeiaTimeStep();
            bool IsInLightCone() const;
            std::pair<unsigned, double> GuessRetardedTime();
            double GetSpaceTimeInterval(const double &aParticleTime, const double &aReceiverTime, const LMCThreeVector &aParticlePosition, const LMCThreeVector &aReceiverPosition ) const;
            double GetStepRoot(const locust::Particle aParticle, double aReceiverTime, LMCThreeVector aReceiverPosition, double aSpaceTimeInterval) const;
            unsigned FindClosestParticle(double tNew) const;
            std::pair<unsigned, double> FindRoot(std::pair<unsigned, double> aRetardedSolution) const;
            void CacheSolution(const int aIndex, const double aRetardedTime);

            locust::Particle fCurrentParticle;
            LMCThreeVector fFieldPosition;
            double fFieldTime;
            double fPreviousFieldTime;
            unsigned fAntennaIndex;
            double fKassiopeiaTimeStep;

            std::vector<LMCThreeVector> fAntennaPositions;
            std::vector<std::pair<unsigned, double> > fCachedSolutions; //Cache the results from previous iteration. [0] is previous index, [1] is corresponding retarded time of previous solution
            std::vector<bool> fHasCachedSolution;
            kl_interface_ptr_t fInterface;

};


} /* namespace locust */

#endif /* LMCLIENARDWIECHERT_HH_ */
