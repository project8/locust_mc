/*
 * LMCChannel.hh
 *
 *  Created on: Mar 1, 2018
 *      Author: nbuzinsky
 */

#ifndef LMCCHANNEL_HH_
#define LMCCHANNEL_HH_
#include <vector>

namespace locust
{
 /*!
 @class Channel
 @author N. Buzinsky
 @brief Used to conveniently group multiple receivers into a channel/ and do phased sum
 @details
 Available configuration options:
 No input parameters
 */
    template <class T>
    class Channel
    {
        public:

            double PhasedSum();
            void AddReceiver(const T receiver){receiverElements.push_back(receiver);}

            T& operator[] (const int &index) {return receiverElements[index];}
            std::size_t size() {return receiverElements.size();}

        private:
            std::vector<T> receiverElements; //Best way to add receivers w/o object slicing?

};


} /* namespace locust */

#endif /* LMCCHANNEL_HH_ */
