/*
 * LMCFieldCalculator.hh
 *
 *  Created on: Oct 4, 2018
 *      Author: pslocum, buzinsky
 */

#ifndef LMCFIELDCALCULATOR_HH_
#define LMCFIELDCALCULATOR_HH_

#include "KSComponentTemplate.h"
#include "KSSpaceInteraction.h"
#include "KToolbox.h"
#include "LMCConst.hh"
#include "KSTrajectory.h"
#include "KSParticle.h"
#include <vector>

namespace locust
{
    /*!
      @class FieldCalculator
      @author P. Slocum
      @brief Base class to compute propagating field amplitudes inside the step modifier LMCCyclotronRadiationExtractor.
      @details
      q Available configuration options:
      No input parameters
      */

    class FieldCalculator 
    {

        public:
            FieldCalculator();

            double GetGroupVelocityTE11(Kassiopeia::KSParticle& aFinalParticle);
            double GetGroupVelocityTM01(Kassiopeia::KSParticle& aFinalParticle);
            double GetGroupVelocityTE01(Kassiopeia::KSParticle& aFinalParticle);
            double GetDampingFactorPhase2(Kassiopeia::KSParticle& anInitialParticle, Kassiopeia::KSParticle& aFinalParticle);
            double GetDampingFactorPhase1(Kassiopeia::KSParticle& anInitialParticle, Kassiopeia::KSParticle& aFinalParticle);
            double GetCouplingFactorTE11(Kassiopeia::KSParticle& aFinalParticle);
            double GetCouplingFactorTM01(Kassiopeia::KSParticle& aFinalParticle);
            double GetCouplingFactorTE01(Kassiopeia::KSParticle& aFinalParticle);
            double GetTM01FieldWithTerminator(Kassiopeia::KSParticle& anInitialParticle, Kassiopeia::KSParticle& aFinalParticle);
            double GetTE11FieldAfterOneBounce(Kassiopeia::KSParticle& anInitialParticle, Kassiopeia::KSParticle& aFinalParticle);
            double GetTE01FieldAfterOneBounce(Kassiopeia::KSParticle& anInitialParticle, Kassiopeia::KSParticle& aFinalParticle);

    };

} /* namespace locust */

#endif /* LMCFIELDCALCULATOR_HH_ */
