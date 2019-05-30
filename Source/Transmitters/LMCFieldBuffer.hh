/*
 * LMCFieldBuffer.hh
 *
 *  Created on: May 19, 2019
 *      Author: pslocum
 */

#ifndef LMCFIELDBUFFER_HH_
#define LMCFIELDBUFFER_HH_

#include <vector>
#include <deque>

namespace locust
{
 /*!
 @class FieldBuffer
 @author P. Slocum
 @brief Class to buffer information about transmitted fields prior to their arrival at the receiver.
 @details
 Available configuration options:
 No input parameters
 */

    class FieldBuffer
    {

        public:
            FieldBuffer();
            virtual ~FieldBuffer();

            std::vector<std::deque<double>> InitializeBuffer(int nchannels, int npatches, int buffersize);

        private:
            std::vector<std::deque<double>> buffer;

    };

} /* namespace locust */

#endif /* LMCPATCHANTENNA_HH_ */
