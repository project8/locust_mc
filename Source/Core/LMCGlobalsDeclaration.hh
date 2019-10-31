/*
 * GlobalsDeclaration.hh
 *
 *  Created on: Sept. 24, 2015
 *      Author: pslocum
 */

#ifndef GLOBALSDECLARATION_HH_
#define GLOBALSDECLARATION_HH_


#include <condition_variable>
#include <vector>
#include <mutex>
#include <deque>
#include "LMCParticle.hh"


extern int Project8Phase; // 1, 2, or 3
extern double CENTER_TO_SHORT;
extern double CENTER_TO_ANTENNA;
extern double t_old;

extern double fKassTimeStep;

extern std::deque<locust::Particle> fParticleHistory;

extern bool fWaitBeforeEvent;
extern bool fWaitAfterEvent;
extern bool fKassEventReady;
extern bool fEventInProgress;
extern bool fRunInProgress;
extern bool fPreEventInProgress;
extern bool fFalseStartKassiopeia;
extern bool fDoneWithSignalGeneration;


extern std::mutex fMutex;
extern std::mutex fKassReadyMutex;
extern std::mutex fMutexDigitizer;


extern std::condition_variable fPreEventCondition;
extern std::condition_variable fPostEventCondition;
extern std::condition_variable fDigitizerCondition;
extern std::condition_variable fKassReadyCondition;



//3 Dimensional arrays: NFDXXXField[ReceiverIndex][Time Series Index][X, Y, Z components]
//It looks bad but is actually optimal: std::arrays give preallocated continuous storage: 
//However, we have indeterminate number of Receiver Points: use std::vector
//std::vector<std::array<std::array<double, 16384>, 3 > >  NFDElectricField;
//std::vector<std::array<std::array<double, 16384>, 3 > >  NFDMagneticField;

#endif /* GLOBALSDECLARATION_HH_ */

