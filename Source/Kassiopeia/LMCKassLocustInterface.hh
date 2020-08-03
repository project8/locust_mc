/*
 * LMCKassLocustInterface.hh
 *
 *  Created on: Aug 21, 2019
 *      Author: N.S. Oblath
 */

#ifndef LOCUST_LMCKASSLOCUSTINTERFACE_HH_
#define LOCUST_LMCKASSLOCUSTINTERFACE_HH_

#include "LMCParticle.hh"

#include "singleton.hh"

#include <condition_variable>
#include <deque>
#include <memory>
#include <mutex>


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
        bool fPreEventInProgress;
        bool fFalseStartKassiopeia; // flag to avoid false start on some Macs.
        bool fDoneWithSignalGeneration;  // do not continue to generate voltages and advance digitizer time.


        std::mutex fMutex;  // pls:  this mutex is used for pre and post event mods.
        std::mutex fKassReadyMutex;
        std::mutex fMutexDigitizer;

        std::condition_variable fPreEventCondition;
        std::condition_variable fPostEventCondition;
        std::condition_variable fDigitizerCondition;
        std::condition_variable fKassReadyCondition;

        // TODO: remove these in some way to avoid the Project 8-specificity
        int fProject8Phase; // 1, 2, or 3, defined with the step modifier instance in the xml file.
        double fCENTER_TO_SHORT;
        double fCENTER_TO_ANTENNA;

    };

    typedef std::shared_ptr< KassLocustInterface > kl_interface_ptr_t;

    class KLInterfaceBootstrapper : public scarab::singleton< KLInterfaceBootstrapper >
    {
        public:
            kl_interface_ptr_t GetInterface() const;

            void SetInterface( kl_interface_ptr_t aInterface );

        protected:
            friend class scarab::singleton< KLInterfaceBootstrapper >;
            friend class scarab::destroyer< KLInterfaceBootstrapper >;

            KLInterfaceBootstrapper();
            ~KLInterfaceBootstrapper();

            kl_interface_ptr_t fInterface;
    };

} /* namespace locust */

#endif /* LOCUST_LMCKASSLOCUSTINTERFACE_HH_ */
