/*
 * LMCTransmitterHardware.cc
 *
 *  Created on: Feb 18, 2020
 *      Author: pslocum
 */

#include "LMCTransmitterHardware.hh"


namespace locust
{
    TransmitterHardware::TransmitterHardware():
    		fDrivePhaseDifference( 0. ),
    		fNAntennas( 1 )
    {

    }
    TransmitterHardware::~TransmitterHardware() {}



    void TransmitterHardware::TxHardwareSayHello()
     {
    	printf("TransmitterHardware says hello\n"); getchar();
     }

    int TransmitterHardware::GetNAntennas()
    {
    	return fNAntennas;
    }

    void TransmitterHardware::SetNAntennas(int aNumber)
    {
    	fNAntennas = aNumber;
    }


} /* namespace locust */

