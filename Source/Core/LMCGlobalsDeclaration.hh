/*
 * GlobalsDeclaration.hh
 *
 *  Created on: Sept. 24, 2015
 *      Author: pslocum
 */

#ifndef GLOBALSDECLARATION_HH_
#define GLOBALSDECLARATION_HH_

#define NPATCHES_PER_STRIP 1
#define PATCH_SPACING 0.0054 // m.  Spacing along Z.
#define PATCH_RADIUS 0.0 // m.  Radius of patch ring.
#define PATCH_RING_OFFSET 0.0 // m.  Offset of entire array in X direction



#include <condition_variable>
#include <vector>
#include <mutex>
#include <deque>
#include "LMCParticle.hh"


extern int Project8Phase; // 1, 2, or 3
extern double CENTER_TO_SHORT;
extern double CENTER_TO_ANTENNA;
extern double t_old;

extern double fDigitizerTimeStep;

extern std::deque<locust::Particle> fParticleHistory;

extern bool fWaitBeforeEvent;
extern bool fWaitAfterEvent;
extern bool fKassEventReady;
extern bool fEventInProgress;
extern bool fRunInProgress;
extern bool fPreEventInProgress;
extern bool fFalseStartKassiopeia;
extern bool fDoneWithSignalGeneration;


extern std::mutex fMutex;  // pls:  this mutex is used for pre and post event mods.
extern std::mutex fKassReadyMutex;  // pls:  this mutex is used for pre and post event mods.
extern std::mutex fMutexDigitizer;  // pls:  not completely sure we need an extra mutex, but it may help clarify.


//  These global condition variables are causing the simulation to hang after writing to the egg file.
extern std::condition_variable fPreEventCondition;
extern std::condition_variable fPostEventCondition;
extern std::condition_variable fDigitizerCondition;
extern std::condition_variable fKassReadyCondition;



//3 Dimensional arrays: NFDXXXField[ReceiverIndex][Time Series Index][X, Y, Z components]
//It looks bad but is actually optimal: std::arrays give preallocated continuous storage: 
//However, we have indeterminate number of Receiver Points: use std::vector
//std::vector<std::array<std::array<double, 16384>, 3 > >  NFDElectricField;  // pls:  placeholder for oversampled signal.
//std::vector<std::array<std::array<double, 16384>, 3 > >  NFDMagneticField;  // pls:  placeholder for oversampled signal.

#endif /* GLOBALSDECLARATION_HH_ */

