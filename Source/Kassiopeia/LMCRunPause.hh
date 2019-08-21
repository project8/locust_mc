/*
 * LMCRunPause.hh
 *
 *  Created on: Jul 31, 2019
 *      Author: N.S. Oblath
 */

#ifndef LOCUST_LMCRUNPAUSE_HH_
#define LOCUST_LMCRUNPAUSE_HH_

#include "KSRunModifier.h"
#include "KSComponentTemplate.h"

#include "LMCKassLocustInterface.hh"

namespace locust
{

    class RunPause :
            public Kassiopeia::KSComponentTemplate< RunPause, Kassiopeia::KSRunModifier >
    {
        public:
            RunPause( kl_interface_ptr_t aInterface );
            RunPause( const RunPause& aOrig );
            virtual ~RunPause();

            RunPause* Clone() const;

        public:

            virtual bool ExecutePreRunModification( Kassiopeia::KSRun& aRun );
            virtual bool ExecutePostRunModification( Kassiopeia::KSRun& aRun );

        protected:
            void WakeAfterEvent(unsigned TotalEvents, unsigned EventsSoFar);
            bool ReceivedEventStartCondition()

        private:
            void InitializeComponent();
            void DeinitializeComponent();

        protected:
            virtual void PullDeupdateComponent();
            virtual void PushDeupdateComponent();

        protected:
            kl_interface_ptr_t fInterface;

    };

} /* namespace locust */

#endif /* LOCUST_LMCRUNPAUSE_HH_ */

