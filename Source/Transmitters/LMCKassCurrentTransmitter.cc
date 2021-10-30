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


    double KassCurrentTransmitter::quadrantOrbitCorrection(double phase, double vx)
    {
    	double phaseCorrection = 0.;
    	if (((phase < 0.)&&(vx < 0.)) || ((phase > 0.)&&(vx > 0.)))
    		phaseCorrection = LMCConst::Pi();

    	return phaseCorrection;
    }


    double KassCurrentTransmitter::quadrantPositionCorrection(double phase, double x)
    {
    	double phaseCorrection = 0.;
    	if (((phase < 0.)&&(x < 0.)) || ((phase > 0.)&&(x < 0.)))
    		phaseCorrection = LMCConst::Pi();

    	return phaseCorrection;
    }


    double KassCurrentTransmitter::calcOrbitPhase(double vx, double vy)
    {
    	double phase = 0.;
    	if (fabs(vy) > 0.)
    		phase = atan(-vx/vy);
    	phase += quadrantOrbitCorrection(phase, vx);
//    	printf("phase is %g\n", phase*180./LMCConst::Pi()); getchar();
    	return phase;
    }

    double KassCurrentTransmitter::calcTheta(double x, double y)
    {
    	double phase = 0.;
    	if (fabs(x) > 0.)
    		phase = atan(y/x);
    	phase += quadrantPositionCorrection(phase, x);
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

        particleXP[0] = pow(tParticle.GetPosition().X()*tParticle.GetPosition().X() + tParticle.GetPosition().Y()*tParticle.GetPosition().Y(), 0.5);
        particleXP[1] = calcTheta(tParticle.GetPosition().X(), tParticle.GetPosition().Y());
        particleXP[2] = tParticle.GetPosition().Z();
        particleXP[3] = tParticle.GetVelocity().X();
        particleXP[4] = tParticle.GetVelocity().Y();
        particleXP[5] = tParticle.GetVelocity().Z();
        particleXP[6] = calcOrbitPhase(tParticle.GetVelocity().X(), tParticle.GetVelocity().Y());
        particleXP[7] = tParticle.GetCyclotronFrequency();

    	return particleXP;
    }


} /* namespace locust */
