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
            fTime( -99. ), 
            fTimeDisplacement( -99.)
    {

    }

    LienardWiechert::~LienardWiechert()
    {

    }


    //Return index of fParticleHistory particle closest to the time we are evaluating
    int LienardWiechert::FindNode(double tNew) const
    {
        std::deque<locust::Particle>::iterator it;

        //Get iterator pointing to particle step closest to tNew
        it = std::upper_bound( fParticleHistory.begin() , fParticleHistory.end() , tNew, [] (const double &a , const locust::Particle &b) { return a < b.GetTime();} );

        int tNodeIndex = it - fParticleHistory.begin();

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

     bool LienardWiechert::IsInLightCone()
    {
        if(fParticleHistory.front().GetTime()<=3.*kassiopeiaTimeStep)
        {
            fParticleHistory.front().Interpolate(0);
            if(GetSpaceTimeInterval(fParticleHistory.front().GetTime(true), tReceiverTime , fParticleHistory.front().GetPosition(true), currentPatch->GetPosition() ) < 0 )
            {
                //printf("Skipping! out of Bounds!: tReceiverTime=%e\n",tReceiverTime);
                continue;
            }
        }
    }

    void LienardWiechert::FindRoot()
    {
        for(int j=0;j<25;++j)
        {

            tRetardedTime = GetStepRoot(tCurrentParticle, tReceiverTime, currentPatch->GetPosition(), tSpaceTimeInterval);
            tCurrentParticle.Interpolate(tRetardedTime);

            //Change the kassiopeia step we expand around if the interpolation time displacement is too large
            if(fabs(tCurrentParticle.GetTime(true) - tCurrentParticle.GetTime(false)) > kassiopeiaTimeStep)
            {
                CurrentIndex=FindNode(tRetardedTime);
                tCurrentParticle=fParticleHistory[CurrentIndex];
                tCurrentParticle.Interpolate(tRetardedTime);
            }

            tSpaceTimeInterval = GetSpaceTimeInterval(tCurrentParticle.GetTime(true), tReceiverTime, tCurrentParticle.GetPosition(true), currentPatch->GetPosition());
            tOldSpaceTimeInterval = tSpaceTimeInterval;
        }
    }

    void LienardWiechert::InitialRetardedTimeGuess()
    {
        if(currentPatch->GetPreviousRetardedIndex() == -99.)
        {
            CurrentIndex=FindNode(tReceiverTime);
            tCurrentParticle = fParticleHistory[CurrentIndex];
            tRetardedTime = tReceiverTime - (tCurrentParticle.GetPosition() - currentPatch->GetPosition() ).Magnitude() / LMCConst::C();
            if(tRetardedTime < 0) 
            {
                tRetardedTime = 0;
                CurrentIndex = 0;
            }
        }
        else
        {
            CurrentIndex = currentPatch->GetPreviousRetardedIndex();
            tRetardedTime = currentPatch->GetPreviousRetardedTime() + tLocustStep;
        }

        CurrentIndex = FindNode(tRetardedTime);
    }

    void LienardWiechert::Cache()
    {
        currentPatch->SetPreviousRetardedIndex(CurrentIndex);
        currentPatch->SetPreviousRetardedTime(tRetardedTime);
    }

} /* namespace locust */
