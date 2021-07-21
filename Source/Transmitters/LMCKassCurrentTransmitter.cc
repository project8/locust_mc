/*
 * LMCKassCurrentTransmitter.cc
 *
 *  Created on: June 10, 2021
 *      Author: pslocum
 */

#include "LMCKassCurrentTransmitter.hh"
#include "logger.hh"
#include <thread>
#include <algorithm>

#include <iostream>
#include <fstream>
#include <math.h>

using std::string;


namespace locust
{
    LOGGER( lmclog, "KassCurrentTransmitter" );

    KassCurrentTransmitter::KassCurrentTransmitter():
    fInterface( KLInterfaceBootstrapper::get_instance()->GetInterface() )
    {
    }

    KassCurrentTransmitter::~KassCurrentTransmitter()
    {
    }

    bool KassCurrentTransmitter::Configure( const scarab::param_node& aParam )
    {
        return true;
    }

    bool KassCurrentTransmitter::IsKassiopeia()
    {
    	return true;
    }


    double KassCurrentTransmitter::quadrantCorrection(double phase, std::vector<double> tKassParticleXP)
    {
    	double phaseCorrection = 0.;
    	if (((phase < 0.)&&(tKassParticleXP[3] < 0.)) || ((phase > 0.)&&(tKassParticleXP[3] > 0.)))
    		phaseCorrection = LMCConst::Pi();

    	return phaseCorrection;
    }

    double KassCurrentTransmitter::calcOrbitPhase(std::vector<double> tKassParticleXP)
    {
    	double phase = 0.;
    	if (fabs(tKassParticleXP[3]) > 0.)
    		phase = atan(-tKassParticleXP[3]/tKassParticleXP[4]);
    	phase += quadrantCorrection(phase, tKassParticleXP);
//    	printf("phase is %g\n", phase*180./LMCConst::Pi()); getchar();
    	return phase;
    }


    //Return index of fParticleHistory particle closest to the time we are evaluating
    int KassCurrentTransmitter::FindNode(double tNew) const
    {
        std::deque<locust::Particle>::iterator it;
        it = std::upper_bound( fInterface->fParticleHistory.begin() , fInterface->fParticleHistory.end() , tNew, [] (const double &a , const locust::Particle &b) { return a < b.GetTime();} );

        int tNodeIndex = it - fInterface->fParticleHistory.begin();

        return tNodeIndex;
    }



    std::vector<double> KassCurrentTransmitter::ExtractParticleXP(double TOld)
    {

        int currentIndex = FindNode(TOld);
        locust::Particle tParticle = fInterface->fParticleHistory[currentIndex];
        tParticle.Interpolate(TOld);

        std::vector<double> particleXP;
    	particleXP.resize(8);

        particleXP[0] = tParticle.GetPosition().X();
        particleXP[1] = tParticle.GetPosition().Y();
        particleXP[2] = tParticle.GetPosition().Z();
        particleXP[3] = tParticle.GetVelocity().X();
        particleXP[4] = tParticle.GetVelocity().Y();
        particleXP[5] = tParticle.GetVelocity().Z();
        particleXP[6] = calcOrbitPhase(particleXP);
        particleXP[7] = tParticle.GetCyclotronFrequency();

    	return particleXP;
    }


} /* namespace locust */
