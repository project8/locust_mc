/*
 * LMCReceiver.cc
 *
 *  Created on: Mar 5, 2018
 *      Author: nbuzinsky
 */

#include "LMCReceiver.hh"

namespace locust
{
    Receiver::Receiver():
		copolarizationDirection(0,0,0),
		crosspolarizationDirection(0,0,0),
		normalDirection(0,0,0),
		centerPosition(0,0,0)
	{}
    Receiver::~Receiver() {}


    LMCThreeVector Receiver::GetPolarizationDirection()
    {
    	return copolarizationDirection;
    }

    void Receiver::SetPolarizationDirection(const LMCThreeVector &copolDirection)
    {
        copolarizationDirection = copolDirection;
    }


    LMCThreeVector Receiver::GetCrossPolarizationDirection()
    {
    	return crosspolarizationDirection;
    }

    void Receiver::SetCrossPolarizationDirection(const LMCThreeVector &crosspolDirection)
    {
        crosspolarizationDirection = crosspolDirection;
    }


    LMCThreeVector Receiver::GetNormalDirection()
    {
    	return normalDirection;
    }

    void Receiver::SetNormalDirection(const LMCThreeVector &normDirection)
    {
        normalDirection = normDirection;
    }

    LMCThreeVector Receiver::GetPosition()
    {
    	return centerPosition;
    }

    void Receiver::SetCenterPosition(const LMCThreeVector &newPosition)
    {
        centerPosition = newPosition;
    }


    void Receiver::RxSayHello()
     {
     	printf("rx says hello\n");
     	getchar();
     }

    double Receiver::GetPatternFactor(LMCThreeVector incidentKVector, Receiver currentElement)
    {
    	return 1.0;
    }





} /* namespace locust */

