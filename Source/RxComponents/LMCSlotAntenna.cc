/*
 * LMCSlotAntenna.cc
 *
 *  Created on: Dec 19, 2019
 *      Author: pslocum
 */

#include "LMCSlotAntenna.hh"


namespace locust
{

    SlotAntenna::SlotAntenna()
    {
    }

    SlotAntenna::~SlotAntenna()
    {
    }

    void SlotAntenna::RxSayHello()
    {
      	printf("slot says hello\n");
      	getchar();
    }

    double SlotAntenna::GetPatternFactor(LMCThreeVector incidentKVector, Receiver currentElement)
    {
    	return 1.0;
    }


} /* namespace locust */
