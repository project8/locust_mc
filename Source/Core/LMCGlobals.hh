/*
 * Globals.hh
 *
 *  Created on: Sept. 24, 2015
 *      Author: pslocum
 */

#ifndef GLOBALS_HH_
#define GLOBALS_HH_

#include <condition_variable>
#include <mutex>

double Z = -99.;
double R = -99.;
double K = -99.;
double vlong = -99.;
double zvelocity = -99.;
double t = 0.;
double t_old = 0.;
double dt = -99.;
double fcyc = -99.;
double de = -99.;
double LarmorPower = -99.;
double GammaZ = -99.;

double testvar = 0.;

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



#endif /* GLOBALS_HH_ */
