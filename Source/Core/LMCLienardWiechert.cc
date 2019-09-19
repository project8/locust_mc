/*
 * LMCLienardWiechert.cc
 *
 *  Created on: Sep 17, 2019
 *      Author: nbuzinsky
 */

#include "LMCLienardWiechert.hh"


namespace locust
{
    LienardWiechert::LienardWiechert() :
            fFieldPosition( 0., 0., 0. ),
            fFieldTime( 0. ),
            fAntennaIndex( 0 ),
            fAntennaPositions(),
            fHasCachedSolution()
    {

    }

    LienardWiechert::~LienardWiechert()
    {

    }

    void LienardWiechert::AddFieldPoint(const LMCThreeVector aFieldPoint)
    {
        fAntennaPositions.push_back(aFieldPoint);
        fPreviousTimes.push_back(std::pair<unsigned, int>(0,0));
        fHasCachedSolution.push_back(false);
    }

    void LienardWiechert::SetFieldEvent(const double aTime, const unsigned aFieldPointIndex)
    {
        fAntennaIndex = aFieldPointIndex;
        fFieldPosition = fAntennaPositions[fAntennaIndex];
        fFieldTime = aTime;
    }

    bool LienardWiechert::SolveFieldSolutions()
    {
        if(!IsInLightCone())
        {
            return false;
        }

        GuessRetardedTime();

        FindRoot();
        CacheSolution(tParticleIndex, tRetardedTime);

        return true;
    }

    LMCThreeVector LienardWiechert::GetElectricField()
    {
        return fCurrentParticle.CalculateElectricField(fFieldPosition);
    }

    LMCThreeVector LienardWiechert::GetMagneticField()
    {
        return fCurrentParticle.CalculateMagneticField(fFieldPosition);
    }


    //Return index of fParticleHistory particle closest to the time we are evaluating
    unsigned LienardWiechert::FindClosestParticle(double tNew) const
    {
        std::deque<locust::Particle>::iterator it;

        //Get iterator pointing to particle step closest to tNew
        it = std::upper_bound( fParticleHistory.begin() , fParticleHistory.end() , tNew, [] (const double &a , const locust::Particle &b) { return a < b.GetTime();} );

        unsigned tNodeIndex = it - fParticleHistory.begin();

        return tNodeIndex;
    }

    double LienardWiechert::GetSpaceTimeInterval(const double &aParticleTime, const double &aReceiverTime, const LMCThreeVector &aParticlePosition, const LMCThreeVector &aReceiverPosition ) const
    {
        return aReceiverTime - aParticleTime - (aReceiverPosition - aParticlePosition).Magnitude() / LMCConst::C();
    }

    double LienardWiechert::GetStepRoot(const locust::Particle aParticle, double aReceiverTime, LMCThreeVector aReceiverPosition, double aSpaceTimeInterval) const
    {
        double tRetardedTime = aParticle.GetTime(true);
        return tRetardedTime + aSpaceTimeInterval;
    }

     bool LienardWiechert::IsInLightCone() const
    {
        if(fParticleHistory.front().GetTime()<=3.*kassiopeiaTimeStep)
        {
            fParticleHistory.front().Interpolate(0);
            if(GetSpaceTimeInterval(fParticleHistory.front().GetTime(true), fAntennaTime , fParticleHistory.front().GetPosition(true), fFieldPosition ) < 0 )
            {
                return false;
            }
        }

        return true;
    }

    std::pair<unsigned, double> LienardWiechert::FindRoot() const
    {
            double tSpaceTimeInterval = GetSpaceTimeInterval(tCurrentParticle.GetTime(true), tFieldTime, tCurrentParticle.GetPosition(true), fFieldPosition);
            double tRetardedTime;
            unsigned tIndex;

        for(int j=0;j<25;++j)
        {
            tRetardedTime = GetStepRoot(tCurrentParticle, fFieldTime, fFieldPosition, tSpaceTimeInterval);
            tCurrentParticle.Interpolate(tRetardedTime);

            //Change the kassiopeia step we expand around if the interpolation time displacement is too large
            if(fabs(tCurrentParticle.GetTime(true) - tCurrentParticle.GetTime(false)) > kassiopeiaTimeStep)
            {
                tIndex = FindNode(tRetardedTime);
                tCurrentParticle = fParticleHistory[tIndex];
                tCurrentParticle.Interpolate(tRetardedTime);
            }

            tSpaceTimeInterval = GetSpaceTimeInterval(tCurrentParticle.GetTime(true), tFieldTime, tCurrentParticle.GetPosition(true), fFieldPosition);
        }

        return std::pair<unsigned, double>(tIndex, tRetardedTime);
    }

    std::pair<unsigned, double> LienardWiechert::GuessRetardedTime()  const
    {
        unsigned tIndex;
        double tRetardedTime;

        if(fHasCachedSolution[fAntennaIndex])
        {
            tIndex = fCachedSolutions[fAntennaIndex].first;
            tRetardedTime = fCachedSolutions[fAntennaIndex].second + tLocustStep;
        }
        else
        {
            tIndex = FindNode(fFieldTime);
            tCurrentParticle = fParticleHistory[tIndex];
            tRetardedTime = tFieldTime - (tCurrentParticle.GetPosition() - fFieldPosition).Magnitude() / LMCConst::C();
            if(tRetardedTime < 0) 
            {
                tRetardedTime = 0;
                tIndex = 0;
            }
        }

        CurrentIndex = FindNode(tRetardedTime);
        return std::pair<unsigned, double>(tIndex, tRetardedTime);
    }

    void LienardWiechert::CacheSolution(const int aCurrentIndex, const double aRetardedTime)
    {
        fPreviousTime[fPatchIndex] = std::pair<int, double>(aCurrentIndex, aRetardedTime);
    }

} /* namespace locust */
