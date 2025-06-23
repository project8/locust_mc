/*
 * LMCRunPause.hh
 *
 *  Created on: Jul 31, 2019
 *      Author: N.S. Oblath
 */

#ifndef LOCUST_LMCRUNPAUSE_HH_
#define LOCUST_LMCRUNPAUSE_HH_

#include "KSRunModifier.h"
#include "KToolbox.h"
#include "KSTermMaxR.h"
#include "KSTermMaxTime.h"
#include "KSTermDeath.h"
#include "KSRootTerminator.h"
#include "KSGeoSpace.h"
#include "KSGeoSurface.h"
#include "KGBoxSpace.hh"
#include "KGCore.hh"
#include "KSRootGenerator.h"
#include "KSGenGeneratorComposite.h"
#include "KSGenValueUniform.h"
#include "KSGenValueAngleSpherical.h"
#include "KSGenValueFix.h"
#include "KSGenEnergyComposite.h"
#include "KSGenEnergyKryptonEvent.h"
#include "KSGenPositionRectangularComposite.h"
#include "KSGenDirectionSphericalComposite.h"
#include "KSGenTimeComposite.h"
#include "LMCKassLocustInterface.hh"
#include "KRandom.h"
#include "LMCDistributionInterface.hh"
#include <stdlib.h>
#include <time.h>


namespace locust
{

    class RunPause :
            public Kassiopeia::KSComponentTemplate< RunPause, Kassiopeia::KSRunModifier >
    {
        public:
            RunPause();
            RunPause( const RunPause& aOrig );
            virtual ~RunPause();

            RunPause* Clone() const;

        public:

            virtual bool ExecutePreRunModification( Kassiopeia::KSRun& aRun );
            virtual bool ExecutePostRunModification( Kassiopeia::KSRun& aRun );

            bool ConfigureByInterface();
            bool Configure( const scarab::param_node& aParam );
            Kassiopeia::KSGeoSpace* GetKSWorldSpace();
            KGeoBag::KGSpace* GetKGWorldSpace();
            bool DeleteLocalKassObjects();
            bool AddWaveguideTerminator( const scarab::param_node& aParam );
            bool AddMaxTimeTerminator( const scarab::param_node& aParam );
            bool AddMaxRTerminator( const scarab::param_node& aParam );
            bool AddGenerator( const scarab::param_node& aParam );
            int GetSeed( const scarab::param_node& aParam );




        protected:
            kl_interface_ptr_t fInterface;

        private:
            katrin::KToolbox& fToolbox;
            Kassiopeia::KSTermMaxTime* fLocustMaxTimeTerminator;
            Kassiopeia::KSTermMaxR* fLocustMaxRTerminator;
            KGeoBag::KGBoxSpace* fBox;
            KGeoBag::KGSpace* fKGSpace;
            Kassiopeia::KSGeoSurface* fSurface;
            Kassiopeia::KSTermDeath* fLocustTermDeath;
            Kassiopeia::KSCommand* fCommand;
            Kassiopeia::KSGeoSpace* fKSSpace;
            Kassiopeia::KSGenGeneratorComposite* fGenerator;
            Kassiopeia::KSGenDirectionSphericalComposite* fGenDirectionComposite;
			Kassiopeia::KSGenValueAngleSpherical* fThetaGenerator;
			Kassiopeia::KSGenValueUniform* fPhiGenerator;
			Kassiopeia::KSGenPositionRectangularComposite* fGenPositionComposite;
			Kassiopeia::KSGenValueUniform* fPositionXGenerator;
			Kassiopeia::KSGenValueUniform* fPositionYGenerator;
			Kassiopeia::KSGenValueUniform* fPositionZGenerator;
			Kassiopeia::KSGenValueUniform* fEnergyUniform;
			Kassiopeia::KSGenEnergyKryptonEvent* fEnergyKrypton;
			Kassiopeia::KSGenEnergyComposite* fGenEnergyComposite;
			Kassiopeia::KSGenCreator* fGenEnergyCreator;
			Kassiopeia::KSGenTimeComposite* fGenTimeComposite;
			Kassiopeia::KSGenValueUniform* fTimeGenerator;
			Kassiopeia::KSGenValueFix* fGenPidComposite;

            std::shared_ptr< BaseDistribution> fTrackLengthDistribution;
            DistributionInterface fDistributionInterface;
            double fMinTrackLengthFraction;
            bool fConfigurationComplete;
            int fEventCounter;
            int fMaxEvents;
            bool fMottScattering;





    };

} /* namespace locust */

#endif /* LOCUST_LMCRUNPAUSE_HH_ */

