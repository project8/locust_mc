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

        double tposX = tParticle.GetPosition(true).X();
        double tposY = tParticle.GetPosition(true).Y();
        double tposZ = tParticle.GetPosition(true).Z();
		double tvX = tParticle.GetVelocity(true).X();
		double tvY = tParticle.GetVelocity(true).Y();
		double tvZ = tParticle.GetVelocity(true).Z();

        std::vector<double> particleXP;
    	particleXP.resize(8);

        particleXP[0] = pow( tposX*tposX + tposY*tposY, 0.5);
        particleXP[1] = calcTheta(tposX, tposY);
        particleXP[2] = tposZ;
        particleXP[3] = tvX;
        particleXP[4] = tvY;
        particleXP[5] = tvZ;
        particleXP[6] = calcOrbitPhase(tvX, tvY);
        particleXP[7] = tParticle.GetCyclotronFrequency();
        particleXP[8] = tParticle.GetLarmorPower();

    	return particleXP;
    }

    void KassCurrentTransmitter::InitializeFieldPoint(LMCThreeVector fieldPoint)
    {
    	Transmitter::InitializeFieldPoint(fieldPoint);
    	fFieldSolver.AddFieldPoint(fieldPoint);
    }

    std::vector<double> KassCurrentTransmitter::SolveCavityKassFields(LMCThreeVector& pointOfInterest, double& tReceiverTime, unsigned& tTotalElementIndex)
    {

	//Currently just solve for fields as if in free space without reflections or mode suppression. Need to update this at some point

        fFieldSolver.SetFieldEvent(tReceiverTime, tTotalElementIndex);
        fFieldSolver.SolveFieldSolutions();

        LMCThreeVector tRadiatedElectricField = fFieldSolver.GetElectricField();
        LMCThreeVector tRadiatedMagneticField = fFieldSolver.GetMagneticField();


        std::vector<double> tSolution;
        tSolution.resize(6);
        tSolution[0] = tRadiatedElectricField.GetX();
        tSolution[1] = tRadiatedElectricField.GetY();
	tSolution[2] = tRadiatedElectricField.GetZ();
	tSolution[3] = tRadiatedMagneticField.GetX();
        tSolution[4] = tRadiatedMagneticField.GetY();
        tSolution[5] = tRadiatedMagneticField.GetZ();
	return tSolution;

    }

} /* namespace locust */
