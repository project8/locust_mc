/*
 * LMCChannel.cc
 *
 *  Created on: Mar 5, 2018
 *      Author: nbuzinsky
 */

#include "LMCReceiver.hh"
#include "LMCChannel.hh"

#include <numeric>

namespace locust
{

    //Use lambda to add all of GetVoltages() from all receiver elements
    template< class T>
    double Channel<T>::PhasedSum()
    {
        return std::accumulate(receiverElements.begin(), receiverElements.end(), 0, 
                               [](double s ,T r) {return s + r->GetVoltage();});
    }

} /* namespace locust */

