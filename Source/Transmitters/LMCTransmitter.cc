/*
 * LMCTransmitter.cc
 *
 *  Created on: Jan 24, 2020
 *      Author: pslocum
 */

#include "LMCTransmitter.hh"


namespace locust
{
    LOGGER( lmclog, "Transmitter" );
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
    
    void Transmitter::SetFieldPoint(int index,LMCThreeVector fieldPoint)
    {
	    fFieldPoints[index]=fieldPoint;
    }
    
    void Transmitter::SetIncidentKVector(int index, LMCThreeVector incidentKVector)
    {
	    fIncidentKVectors[index]=incidentKVector;
    }
    
    void Transmitter::SetPropagationPhaseDelay(int index,double phaseDelay)
    {
	    fPropagationPhaseDelays[index]=phaseDelay;
    }
	   
    void Transmitter::AddPropagationPhaseDelay(LMCThreeVector fieldPoint)
    {
	    fPropagationPhaseDelays.push_back(0.0);
    }

    void Transmitter::AddPropagationPhaseDelay(double phaseDelay)
    {
	    fPropagationPhaseDelays.push_back(phaseDelay);
    }

    void Transmitter::InitializeFieldPoint(LMCThreeVector fieldPoint)
    {
	    AddFieldPoint(fieldPoint);
	    AddIncidentKVector(fieldPoint);
	    AddPropagationPhaseDelay(fieldPoint);
    }

    double Transmitter::GetMeanofFieldPoints(int axis)
    {
	int nPoints=GetNPoints();
	double meanvalue=0.0;
	for(int i=0;i<nPoints;++i)
	{
	     switch (axis)
	     {
		 case 0:
		     meanvalue+=GetFieldPoint(i).GetX();
		     break;
		 case 1:
		     meanvalue+=GetFieldPoint(i).GetY();
		     break;
		 case 2:
		     meanvalue+=GetFieldPoint(i).GetZ();
		     break;
		 default:
		     LWARN(lmclog,"Wrong axis chosen for calculating the mean of field points");
		     LWARN(lmclog,"Please choose between 0,1 or 2 for X,Y, or Z respectively");
		     break;
	     }
	}	
	meanvalue=meanvalue/nPoints;
	return meanvalue;
    }
    
    LMCThreeVector Transmitter::GetMeanofFieldPoints()
    {
	LMCThreeVector meanThreeVector(GetMeanofFieldPoints(0),GetMeanofFieldPoints(1),GetMeanofFieldPoints(2));
	return meanThreeVector;
    }

    LMCThreeVector Transmitter::GetFieldPoint(int index)
    {
	    return fFieldPoints.at(index); 
    }

    LMCThreeVector Transmitter::GetIncidentKVector(int index)
    {
	    return fIncidentKVectors.at(index); 
    }

    double Transmitter::GetPropagationPhaseDelay(int index)
    {
	    return fPropagationPhaseDelays.at(index);
    }

    double Transmitter::GetNPoints()
    {
	    return fFieldPoints.size();
    }
} /* namespace locust */

