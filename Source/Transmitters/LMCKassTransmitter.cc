/*
 * LMCKassTransmitter.cc
 *
 *  Created on: May 4, 2019
 *      Author: pslocum
 */

#include "LMCKassTransmitter.hh"
#include "logger.hh"
#include <thread>
#include <algorithm>

#include <iostream>
#include <fstream>
#include <math.h>

using std::string;


namespace locust
{
    LOGGER( lmclog, "KassTransmitter" );

    KassTransmitter::KassTransmitter()
    {
    }

    KassTransmitter::~KassTransmitter()
    {
    }

    bool KassTransmitter::Configure( const scarab::param_node& aParam )
    {

        return true;
    }

    bool KassTransmitter::IsKassiopeia()
    {
    	return true;
    }



    double* KassTransmitter::GetEFieldCoPol(Receiver* currentElement, int z_index, double elementSpacing, int nElementsPerStrip, double dt)
    {
        double* tSolution = new double[2];
	    return tSolution;

    }



} /* namespace locust */
