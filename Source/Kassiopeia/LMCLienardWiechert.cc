/*
 * LMCLienardWiechert.cc
 *
 *  Created on: Sep 17, 2019
 *      Author: nbuzinsky
 */

#include <deque>

#include "LMCLienardWiechert.hh"

namespace locust
{
    LienardWiechert::LienardWiechert() :
    		fFieldPosition( 0., 0., 0. ),
            fFieldTime( 0. ),
            fPreviousFieldTime( 0. ),
            fAntennaIndex( 0 ),
            fKassiopeiaTimeStep( 0. ),
            fAntennaPositions(),
            fCurrentParticle(),
            fHasCachedSolution(),
            fInterface( KLInterfaceBootstrapper::get_instance()->GetInterface() )

    {
    }

    LienardWiechert::~LienardWiechert()
    {

    }

    void LienardWiechert::SetParticleHistory( std::deque<Particle> aParticleHistory )
    {
    	fInterface->fParticleHistory = aParticleHistory;
    }

    void LienardWiechert::AddFieldPoint(const LMCThreeVector aFieldPoint)
    {
        fAntennaPositions.push_back(aFieldPoint);
        fCachedSolutions.push_back(std::pair<unsigned, int>(0,0));
        fHasCachedSolution.push_back(false);
    }

    void LienardWiechert::SetFieldEvent(const double aTime, const unsigned aFieldPointIndex)
    {
        fAntennaIndex = aFieldPointIndex;
        fFieldPosition = fAntennaPositions[fAntennaIndex];
        fPreviousFieldTime = fFieldTime;
        fFieldTime = aTime;
    }

    void LienardWiechert::SolveFieldSolutions()
    {
        if(!fKassiopeiaTimeStep)
            SetKassiopeiaTimeStep();

        if(!IsInLightCone())
        {
            return;
        }

        std::pair<unsigned, double> tRetardedSolution;
        tRetardedSolution = GuessRetardedTime();

        tRetardedSolution = FindRoot(tRetardedSolution);
        unsigned tIndex = tRetardedSolution.first;
        double tRetardedTime = tRetardedSolution.second;

        fCurrentParticle = fInterface->fParticleHistory[tIndex];
        fCurrentParticle.Interpolate(tRetardedTime);

        CacheSolution(tIndex, tRetardedTime);

    }


    LMCThreeVector LienardWiechert::GetElectricField() const
    {
        if(IsInLightCone())
            return fCurrentParticle.CalculateElectricField(fFieldPosition);
        else
            return LMCThreeVector(0,0,0);
    }

    LMCThreeVector LienardWiechert::GetMagneticField() const
    {
        if(IsInLightCone())
            return fCurrentParticle.CalculateMagneticField(fFieldPosition);
        else
            return LMCThreeVector(0,0,0);
    }

    locust::Particle LienardWiechert::GetRetardedParticle() const
    {
        return fCurrentParticle;
    }

    //Return index of fParticleHistory particle closest to the time we are evaluating
    unsigned LienardWiechert::FindClosestParticle(double tNew) const
    {
        std::deque<locust::Particle>::iterator it;

        //Get iterator pointing to particle step closest to tNew
        it = std::upper_bound( fInterface->fParticleHistory.begin() , fInterface->fParticleHistory.end() , tNew, [] (const double &a , const locust::Particle &b) { return a < b.GetTime();} );

        unsigned tNodeIndex = it - fInterface->fParticleHistory.begin();

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

    void LienardWiechert::SetKassiopeiaTimeStep()
    {
        fKassiopeiaTimeStep = fabs(fInterface->fParticleHistory[0].GetTime() - fInterface->fParticleHistory[1].GetTime());
    }

     bool LienardWiechert::IsInLightCone() const
    {

        if(fInterface->fParticleHistory.front().GetTime()<=3.*fKassiopeiaTimeStep)
        {
            fInterface->fParticleHistory.front().Interpolate(0);
            if(GetSpaceTimeInterval(fInterface->fParticleHistory.front().GetTime(true), fFieldTime , fInterface->fParticleHistory.front().GetPosition(true), fFieldPosition ) < 0 )
            {
                return false;
            }
        }

        return true;
    }

    std::pair<unsigned, double> LienardWiechert::FindRoot(std::pair<unsigned, double> aRetardedSolution) const
    {
        unsigned tIndex = aRetardedSolution.first;
        double tRetardedTime = aRetardedSolution.second;

        locust::Particle tCurrentParticle = fInterface->fParticleHistory[tIndex];
        double tSpaceTimeInterval = GetSpaceTimeInterval(tCurrentParticle.GetTime(true), fFieldTime, tCurrentParticle.GetPosition(true), fFieldPosition);

        for(int j=0;j<25;++j)
        {
            tRetardedTime = GetStepRoot(tCurrentParticle, fFieldTime, fFieldPosition, tSpaceTimeInterval);
            tCurrentParticle.Interpolate(tRetardedTime);

            //Change the kassiopeia step we expand around if the interpolation time displacement is too large
            if(fabs(tCurrentParticle.GetTime(true) - tCurrentParticle.GetTime(false)) > fKassiopeiaTimeStep)
            {
                tIndex = FindClosestParticle(tRetardedTime);
                tCurrentParticle = fInterface->fParticleHistory[tIndex];
                tCurrentParticle.Interpolate(tRetardedTime);
            }

            tSpaceTimeInterval = GetSpaceTimeInterval(tCurrentParticle.GetTime(true), fFieldTime, tCurrentParticle.GetPosition(true), fFieldPosition);
        }

        return std::pair<unsigned, double>(tIndex, tRetardedTime);
    }

    std::pair<unsigned, double> LienardWiechert::GuessRetardedTime()
    {
        unsigned tIndex;
        double tRetardedTime;


        if(fHasCachedSolution[fAntennaIndex])
        {
            tIndex = fCachedSolutions[fAntennaIndex].first;
            tRetardedTime = fCachedSolutions[fAntennaIndex].second + (fFieldTime - fPreviousFieldTime);
        }
        else
        {
            tIndex = FindClosestParticle(fFieldTime);
            locust::Particle tCurrentParticle = fInterface->fParticleHistory[tIndex];
            tRetardedTime = fFieldTime - (tCurrentParticle.GetPosition() - fFieldPosition).Magnitude() / LMCConst::C();
            if(tRetardedTime < 0) 
            {
                tRetardedTime = 0;
                tIndex = 0;
            }
            fHasCachedSolution[fAntennaIndex] = true;
        }

        tIndex = FindClosestParticle(tRetardedTime);
        return std::pair<unsigned, double>(tIndex, tRetardedTime);
    }

    void LienardWiechert::CacheSolution(const int aIndex, const double aRetardedTime)
    {
        fCachedSolutions[fAntennaIndex] = std::pair<int, double>(aIndex, aRetardedTime);
    }

} /* namespace locust */
