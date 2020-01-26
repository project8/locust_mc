/*
 * LMCTransmitter.hh
 *
 *  Created on: Jan. 24, 2020
 *      Author: pslocum
 */

#ifndef LMCTRANSMITTER_HH_
#define LMCTRANSMITTER_HH_


namespace locust
{
 /*!
 @class Transmitter
 @author P. Slocum
 @brief Base class to characterize Transmitter selection
 @details
 Available configuration options:
 No input parameters
 */


    class Transmitter
    {

        public:
            Transmitter();
            Transmitter(double testvar);
            virtual ~Transmitter();
            virtual void TxSayHello();


        private:


};


} /* namespace locust */

#endif /* LMCCHANNEL_HH_ */
