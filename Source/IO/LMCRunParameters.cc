/*
 * LMCRunParameters.cc
 *
 * This class has access to both locust and ROOT libraries.  The syntax is
 * consistent with KTROOTData.cc and the instructions in
 * https://root.cern.ch/root/Using.html .   It is also mentioned in LMCEventLinkDef.hh .
 *
 *  Created on: Jul 7, 2020
 *      Author: pslocum
 */


#include "LMCRunParameters.hh"
#include <iostream>

ClassImp(locust::RunParameters);

namespace locust
{
    RunParameters::RunParameters() {}
    RunParameters::~RunParameters() {}

}
