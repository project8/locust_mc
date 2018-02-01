/*
 * GlobalsDefinition.hh
 *
 *  Created on: Sept. 24, 2015
 *      Author: pslocum
 */

#ifndef GLOBALSDEFINITION_HH_
#define GLOBALSDEFINITION_HH_

double t_old = -99.;
double t_poststep = -99.;

double fDigitizerTimeStep = 5e-10; //Time step for sampling

//running deque for saving previous few ns of particle history 
//in order to caluclate retarded fields
std::deque<locust::Particle> fParticleHistory;
std::deque<locust::Particle> fNewParticleHistory;

bool fWaitBeforeEvent = true;
bool fWaitAfterEvent = true;
bool fKassEventReady = false;
bool fEventInProgress = false;
bool fRunInProgress = false;
bool fPreEventInProgress = false;
bool fFalseStartKassiopeia = true; // flag to avoid false start on some Macs.

std::mutex fMutex;  // pls:  this mutex is used for pre and post event mods.
std::mutex fKassReadyMutex;  
std::mutex fMutexDigitizer;  // pls:  not completely sure we need an extra mutex, but it may help clarify.

std::condition_variable fPreEventCondition;
std::condition_variable fPostEventCondition;
std::condition_variable fDigitizerCondition;
std::condition_variable fKassReadyCondition;

#endif /* GLOBALSDEFINITION_HH_ */

