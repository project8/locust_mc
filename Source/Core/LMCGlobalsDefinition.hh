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

double phi_shortTE11 = 0.;
double* phi_shortTM01 = new double[200];
double* phi_polarizerTM01 = new double[200];
double testvar = -99.;
double fDigitizerTimeStep = 5e-10; //Time step for sampling
//double fEventModTimeStep = 5e-10; //Exact difference in time for fNewParticle History. 

int fDecimationFactor=10.;

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

//If true Phase II (kass build) if false Phase III (freefield gnerator)
bool fPhaseIISimulation = true;


std::mutex fMutex;  // pls:  this mutex is used for pre and post event mods.
std::mutex fKassReadyMutex;  // pls:  this mutex is used for pre and post event mods.
std::mutex fMutexDigitizer;  // pls:  not completely sure we need an extra mutex, but it may help clarify.

std::condition_variable fPreEventCondition;
std::condition_variable fPostEventCondition;
std::condition_variable fDigitizerCondition;
std::condition_variable fKassReadyCondition;

//double* aLongSignal = new double[41943040];  // pls:  placeholder for oversampled signal.


#endif /* GLOBALSDEFINITION_HH_ */

