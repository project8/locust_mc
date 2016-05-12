/*
 * GlobalsDefinition.hh
 *
 *  Created on: Sept. 24, 2015
 *      Author: pslocum
 */

#ifndef GLOBALSDEFINITION_HH_
#define GLOBALSDEFINITION_HH_

double Z = -99.;
double X = -99.;
double t_poststep = -99.;
double t_old = -99.;
double LarmorPower = -99.;
double xvelocity = -99.;
double yvelocity = -99.;
double zvelocity = -99.;
double fcyc = -99.;
double GammaZ = -99.;
double testvar = -99.;


bool fWaitBeforeEvent = true;
bool fWaitAfterEvent = true;
bool fKassEventReady = false;
bool fEventInProgress = false;
bool fRunInProgress = false;
bool fPreEventInProgress = false;



std::mutex fMutex;  // pls:  this mutex is used for pre and post event mods.
std::mutex fMutexDigitizer;  // pls:  not completely sure we need an extra mutex, but it may help clarify.



//  These global condition variables are causing the simulation to hang after writing to the egg file.
std::condition_variable fPreEventCondition;
std::condition_variable fPostEventCondition;
std::condition_variable fDigitizerCondition;



#endif /* GLOBALSDEFINITION_HH_ */
