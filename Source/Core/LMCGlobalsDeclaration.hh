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
#include <vector>
#include <mutex>
#include <deque>
#include "LMCParticleSlim.hh"


extern double t_poststep;
extern double t_old;
extern double testvar;
extern double EventModTimeStep;

extern int fDecimationFactor;

extern std::deque<locust::ParticleSlim> fParticleHistory;
extern std::deque<locust::ParticleSlim> fNewParticleHistory;

extern bool fWaitBeforeEvent;
extern bool fWaitAfterEvent;
extern bool fKassEventReady;
extern bool fEventInProgress;
extern bool fRunInProgress;
extern bool fPreEventInProgress;
extern bool fFalseStartKassiopeia;


extern std::mutex fMutex;  // pls:  this mutex is used for pre and post event mods.
extern std::mutex fKassReadyMutex;  // pls:  this mutex is used for pre and post event mods.
extern std::mutex fMutexDigitizer;  // pls:  not completely sure we need an extra mutex, but it may help clarify.


//  These global condition variables are causing the simulation to hang after writing to the egg file.
extern std::condition_variable fPreEventCondition;
extern std::condition_variable fPostEventCondition;
extern std::condition_variable fDigitizerCondition;
extern std::condition_variable fKassReadyCondition;


extern double* aLongSignal;  // pls:  placeholder for oversampled signal.

//3 Dimensional arrays: NFDXXXField[ReceiverIndex][Time Series Index][X, Y, Z components]
//It looks bad but is actually optimal: std::arrays give preallocated continuous storage: 
//However, we have indeterminate number of Receiver Points: use std::vector
//std::vector<std::array<std::array<double, 16384>, 3 > >  NFDElectricField;  // pls:  placeholder for oversampled signal.
//std::vector<std::array<std::array<double, 16384>, 3 > >  NFDMagneticField;  // pls:  placeholder for oversampled signal.

#endif /* GLOBALSDECLARATION_HH_ */

