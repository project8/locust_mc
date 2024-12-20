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
    fOrbitPhase( 0. ),
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

        if (phase < 0.)
            phaseCorrection = LMCConst::Pi() + (x > 0.) * LMCConst::Pi();
       	else
       	    phaseCorrection = (x < 0.) * LMCConst::Pi();

        return phaseCorrection;
    }


    double KassCurrentTransmitter::calcOrbitPhase(double vx, double vy)
    {
    	double phase = 0.;
        if ((fabs(vy) > 0.))
    	{
    		phase = atan(-vx/vy);
    	}

    	phase += quadrantOrbitCorrection(phase, vx);
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


    double KassCurrentTransmitter::GetGuidingCenterVy()
    {
    	int dSize = fInterface->fParticleHistory.size();
    	if (dSize > 1)
    	{
    		locust::Particle tParticleSecondToLast = fInterface->fParticleHistory.at(dSize-2);
    		locust::Particle tParticleLast = fInterface->fParticleHistory.back();
    		double tDisplacement = tParticleLast.GetGuidingCenterPosition().Y() - tParticleSecondToLast.GetGuidingCenterPosition().Y();
    		double tTimeInterval = tParticleLast.GetTime() - tParticleSecondToLast.GetTime();
    		double tGuidingCenterVy = tDisplacement / tTimeInterval;
    		return tGuidingCenterVy;
    	}
    	else
    	{
    		return 0.;
    	}
    }


    //Return index of fParticleHistory particle closest to the time we are evaluating
    int KassCurrentTransmitter::FindNode(double tNew) const
    {
        std::deque<locust::Particle>::iterator it;
        it = std::upper_bound( fInterface->fParticleHistory.begin() , fInterface->fParticleHistory.end() , tNew, [] (const double &a , const locust::Particle &b) { return a < b.GetTime();} );

        int tNodeIndex = it - fInterface->fParticleHistory.begin();

        return tNodeIndex;
    }



    std::vector<double> KassCurrentTransmitter::ExtractParticleXP(double TOld, bool Interpolate, double dt, bool bWaveguide)
    {

    	locust::Particle tParticle;
    	std::vector<double> particleXP;

    	int pSize = fInterface->fParticleHistory.size();
    	if ( pSize > 0 )
    	{
    	    if (!Interpolate)
    	    {
    		    tParticle = fInterface->fParticleHistory.back();
    	    }
    	    else
    	    {
    		    int currentIndex = FindNode(TOld);
    		    tParticle = fInterface->fParticleHistory[currentIndex];
    		    tParticle.Interpolate(TOld);
    	    }

            double tposX = tParticle.GetPosition(true).X();
            double tposY = tParticle.GetPosition(true).Y();
            double tposZ = tParticle.GetPosition(true).Z();
            double tvX = tParticle.GetVelocity(true).X();
            double tvY = tParticle.GetVelocity(true).Y();
            double tvZ = tParticle.GetVelocity(true).Z();

        	if (!(std::isfinite(tvX)))  // Check for discontinuity at scattering interaction.
        	{
        		tvX = 0.;
        		tvY = 0.;
        		tvZ = 0.;
        	}


//            fOrbitPhase += dt*tParticle.GetCyclotronFrequency();  // accumulate phase semi-analytically.
    	    fOrbitPhase = calcOrbitPhase(tvX, tvY); // or use Kass kinematics

            particleXP.push_back(pow( tposX*tposX + tposY*tposY, 0.5));
            particleXP.push_back(calcTheta(tposX, tposY));
            particleXP.push_back(tposZ);
            particleXP.push_back(tvX);
            particleXP.push_back(tvY);
            particleXP.push_back(tvZ);
            particleXP.push_back(fOrbitPhase);
            particleXP.push_back(tParticle.GetCyclotronFrequency());
            particleXP.push_back(tParticle.GetLarmorPower());
            particleXP.push_back(tParticle.GetTime());


            if (bWaveguide)  // rotate from Kass to Locust coordinate system?
            {
            	/*
                particleXP[0] = pow( tposZ*tposZ + tposX*tposX, 0.5);
                particleXP[1] = calcTheta(tposZ, tposX);
                particleXP[2] = tposY;  // z->y
                particleXP[3] = tvY;    // x->z
                particleXP[4] = tvX;    // y->x
                particleXP[5] = GetGuidingCenterVy();    // z->y, waveguide axis.
                */
            }

    	}

    	return particleXP;
    }




} /* namespace locust */
