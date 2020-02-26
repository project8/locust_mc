/*
 * LMCTransmitter.cc
 *
 *  Created on: Jan 24, 2020
 *      Author: pslocum
 */

#include "LMCTransmitter.hh"


namespace locust
{
    Transmitter::Transmitter() {}
    Transmitter::~Transmitter() {}



    void Transmitter::TxSayHello()
     {

     }

    void Transmitter::AddFieldPoint(LMCThreeVector fieldPoint)
    {
	    fFieldPoints.push_back(fieldPoint);
    }

    void Transmitter::AddIncidentKVector(LMCThreeVector incidentKVector)
    {
	    fIncidentKVectors.push_back(incidentKVector);
    }

    void Transmitter::SetIncidentKVector(int index, LMCThreeVector incidentKVector)
    {
	    fIncidentKVectors[index]=incidentKVector;
    }

    void Transmitter::InitializeFieldPoint(LMCThreeVector fieldPoint)
    {
	    AddFieldPoint(fieldPoint);
	    AddIncidentKVector(fieldPoint);
    }

    LMCThreeVector Transmitter::GetFieldPoint(int index)
    {
	    return fFieldPoints.at(index); 
    }

    LMCThreeVector Transmitter::GetIncidentKVector(int index)
    {
	    return fIncidentKVectors.at(index); 
    }

} /* namespace locust */

