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

    KassCurrentTransmitter::KassCurrentTransmitter()
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

    std::vector<double> KassCurrentTransmitter::ExtractParticleXP()
    {
        locust::Particle tParticle = fFieldSolver.GetInstantaneousParticle();
    	std::vector<double> particleXP;
    	particleXP.resize(6);

        particleXP[0] = tParticle.GetPosition().X();
        particleXP[1] = tParticle.GetPosition().Y();
        particleXP[2] = tParticle.GetPosition().Z();
        particleXP[3] = tParticle.GetVelocity().X();
        particleXP[4] = tParticle.GetVelocity().Y();
        particleXP[5] = tParticle.GetVelocity().Z();

    	return particleXP;
    }


} /* namespace locust */
