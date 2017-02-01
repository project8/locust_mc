/*
 * GlobalsDeclaration.hh
 *
 *  Created on: Sept. 24, 2015
 *      Author: pslocum
 */

#ifndef GLOBALSDECLARATION_HH_
#define GLOBALSDECLARATION_HH_
#define PI 3.1415926

#include <condition_variable>
#include <mutex>
#include <deque>
#include "LMCParticleSlim.hh"


extern double t_poststep;
extern double t_old;
//extern double LarmorPower;
//extern double X;
//extern double Y;
//extern double Z;
//extern double xVelocity;
//extern double yVelocity;
//extern double zVelocity;
//extern double xMagneticField;
//extern double yMagneticField;
//extern double zMagneticField;
//extern double mparticle;
//extern double qparticle;
//extern double fcyc;
//extern double GammaZ;
extern double testvar;
extern double EventModTimeStep;

extern std::deque<locust::ParticleSlim> fParticleHistory;


extern bool fWaitBeforeEvent;
extern bool fWaitAfterEvent;
extern bool fKassEventReady;
extern bool fEventInProgress;
extern bool fRunInProgress;
extern bool fPreEventInProgress;



extern std::mutex fMutex;  // pls:  this mutex is used for pre and post event mods.
extern std::mutex fMutexDigitizer;  // pls:  not completely sure we need an extra mutex, but it may help clarify.


//  These global condition variables are causing the simulation to hang after writing to the egg file.
extern std::condition_variable fPreEventCondition;
extern std::condition_variable fPostEventCondition;
extern std::condition_variable fDigitizerCondition;
extern std::condition_variable fKassReadyCondition;


extern double* aLongSignal;  // pls:  placeholder for oversampled signal.

#endif /* GLOBALSDECLARATION_HH_ */
