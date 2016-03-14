/*
 * LMCCyclotronRadiationExtractor.hh
 *
 *  Created on: Mar 13, 2016
 *      Author: nsoblath
 */

#ifndef LOCUST_LMCCYCLOTRONRADIATIONEXTRACTOR_HH_
#define LOCUST_LMCCYCLOTRONRADIATIONEXTRACTOR_HH_

#include "KField.h"
#include "KSStepModifier.h"
#include "KSComponentTemplate.h"

namespace locust
{

    class CyclotronRadiationExtractor :
            public Kassiopeia::KSComponentTemplate< CyclotronRadiationExtractor, Kassiopeia::KSStepModifier >
    {
        public:
            CyclotronRadiationExtractor();
            CyclotronRadiationExtractor( const CyclotronRadiationExtractor& aOrig );
            virtual ~CyclotronRadiationExtractor();

            CyclotronRadiationExtractor* Clone() const;

        public:
            bool ExecutePreStepModification( Kassiopeia::KSParticle& anInitialParticle, Kassiopeia::KSParticleQueue& aQueue );
            bool ExecutePostStepModifcation( Kassiopeia::KSParticle& anInitialParticle, Kassiopeia::KSParticle& aFinalParticle, Kassiopeia::KSParticleQueue& aQueue );


        private:
            void InitializeComponent();
            void DeinitializeComponent();

        protected:
            virtual void PullDeupdateComponent();
            virtual void PushDeupdateComponent();

    };

} /* namespace locust */

#endif /* LOCUST_LMCCYCLOTRONRADIATIONEXTRACTOR_HH_ */
