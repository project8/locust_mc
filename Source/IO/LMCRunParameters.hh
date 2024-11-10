/*
 * LMCRunParameters.hh
 *
 * This class has access to both locust and ROOT libraries.  The syntax is
 * consistent with KTROOTData.hh and the instructions in
 * https://root.cern.ch/root/Using.html .  It is also mentioned in LMCEventLinkDef.hh .
 *  Created on: Jul 7, 2020
 *      Author: pslocum
 */



#ifndef LMCRUNPARAMETERS_HH_
#define LMCRUNPARAMETERS_HH_

#include "TObject.h"
#include "LMCRunParameters.hh"
#include "time.h"
#include <sys/time.h>

namespace locust
{

    class RunParameters : public TObject
    {

        public:
            RunParameters();
            virtual ~RunParameters();

            double fNoise;
            double fLOfrequency;
            double fSamplingRateMHz;
            double fDecimationFactor;
            std::string fDataType;
            std::string fSimulationType;
            std::string fSimulationSubType;
            long int fRunID;

            ClassDef(RunParameters,1)  // Root syntax.

    };

}
#endif /* LMCRUNPARAMETERS_HH_ */
