/*
 * LMCEventHold.hh
 *
 *  Created on: Mar 13, 2016
 *      Author: nsoblath
 */

#ifndef LOCUST_LMCEVENTHOLD_HH_
#define LOCUST_LMCEVENTHOLD_HH_

#include "KSEventModifier.h"
#include "KSComponentTemplate.h"


#include "LMCKassLocustInterface.hh"

#ifdef ROOT_FOUND
    #include "LMCRootTreeWriter.hh"
#endif


namespace locust
{

    class EventHold :
            public Kassiopeia::KSComponentTemplate< EventHold, Kassiopeia::KSEventModifier >
    {
        public:
            EventHold();
            EventHold( const EventHold& aOrig );
            bool OpenEvent();
            bool OpenFile();
            bool WriteEvent();
            bool WriteJsonFile();
            virtual ~EventHold();

            EventHold* Clone() const;
            bool ConfigureByInterface();
            bool Configure( const scarab::param_node& aParam );
            std::string fTruthOutputFilename;
            std::string fJsonFileName;
            bool fAccumulateTruthInfo;


        public:

            virtual bool ExecutePreEventModification(Kassiopeia::KSEvent &anEvent);
            virtual bool ExecutePostEventModification(Kassiopeia::KSEvent &anEvent);

        protected:
            kl_interface_ptr_t fInterface;

        private:
            bool fConfigurationComplete;
            long int fEventSeed;
            double fConfiguredEMin;
            double fConfiguredPitchMin;
            double fConfiguredXMin;
            int fEventCounter;

    };

} /* namespace locust */

#endif /* LOCUST_LMCEVENTHOLD_HH_ */

