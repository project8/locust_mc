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





} /* namespace locust */

