/*
 * LMCOpenWG.hh
 *
 *  Created on: Apr 6, 2022
 *      Author: P. T. Surukuchi
 */

#ifndef LMCOpenWG_HH_
#define LMCOpenWG_HH_

#include "LMCTransmitterHardware.hh"


namespace locust
{
 /*!
 @class OpenWG
 @author P. T. Surukuchi
 @brief Derived class describing the orientation and angular dependence of a open waveguide antenna transmitter.
 @details
 Available configuration options:
 No input parameters
 */
    class OpenWG: public TransmitterHardware
    {

        public:
            OpenWG();
            virtual ~OpenWG();

            virtual bool Configure( const scarab::param_node& aNode );

            virtual void TxHardwareSayHello();

            virtual double GetPatternFactor(LMCThreeVector pointOfInterest, int antennaNumber);

        private:

            LMCThreeVector fMomentVector;


    };


} /* namespace locust */

#endif /* LMCOpenWG_HH_ */
