/*
 * Globals.hh
 *
 *  Created on: Sept. 24, 2015
 *      Author: pslocum
 */

#ifndef GLOBALS_HH_
#define GLOBALS_HH_

#include <pthread.h>

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
pthread_mutex_t mymutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t tick;


#endif /* GLOBALS_HH_ */
