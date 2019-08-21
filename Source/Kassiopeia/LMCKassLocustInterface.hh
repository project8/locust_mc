/*
 * LMCKassLocustInterface.hh
 *
 *  Created on: Aug 21, 2019
 *      Author: N.S. Oblath
 */

#ifndef LOCUST_LMCKASSLOCUSTINTERFACE_HH_
#define LOCUST_LMCKASSLOCUSTINTERFACE_HH_

#include <condition_variable>
#include <deque>
#include <memory>
#include <mutex>
#include "LMCParticle.hh"


namespace locust
{

    struct KassLocustInterface
    {
        KassLocustInterface();

        double fTOld;
        double fKassTimeStep; //Time step for sampling

        //running deque for saving previous few ns of particle history
        //in order to caluclate retarded fields
        std::deque<Particle> fParticleHistory;

        bool fWaitBeforeEvent;
        bool fWaitAfterEvent;
        bool fKassEventReady;
        bool fEventInProgress;
        bool fRunInProgress;
        bool fPreEventInProgress;
        bool fFalseStartKassiopeia; // flag to avoid false start on some Macs.
        bool fDoneWithSignalGeneration;  // do not continue to generate voltages and advance digitizer time.


        std::mutex fMutex;  // pls:  this mutex is used for pre and post event mods.
        std::mutex fKassReadyMutex;
        std::mutex fMutexDigitizer;  // pls:  not completely sure we need an extra mutex, but it may help clarify.

        std::condition_variable fPreEventCondition;
        std::condition_variable fPostEventCondition;
        std::condition_variable fDigitizerCondition;
        std::condition_variable fKassReadyCondition;
    };

    typedef std::shared_ptr< KassLocustInterface > kl_interface_ptr_t;

} /* namespace locust */

#endif /* LOCUST_LMCKASSLOCUSTINTERFACE_HH_ */
