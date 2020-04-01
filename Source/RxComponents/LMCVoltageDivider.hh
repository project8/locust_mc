/*
 * LMCVoltageDivider.hh
 *
 *  Created on: Feb. 28, 2020
 *      Author: pslocum
 */

#ifndef LMCVOLTAGEDIVIDER_HH_
#define LMCVOLTAGEDIVIDER_HH_

#include "LMCPowerCombiner.hh"
#include "param.hh"
#include "logger.hh"
#include <iostream>

namespace locust
{
 /*!
 @class VoltageDivider
 @author P. Slocum
 @brief Derived class describing an energy-conserving voltage divider.
 @details
 Available configuration options:
 No input parameters
 */


    class VoltageDivider: public PowerCombiner
    {

        public:
            VoltageDivider();
            virtual ~VoltageDivider();
            virtual bool Configure( const scarab::param_node& aNode );


            std::vector<double> GetResistances(double RJunction, double R0, double RGround, int NPAIRS);
            std::vector<double> GetPartialGains(double RJunction, double R0, double RGround, int NPAIRS);
            double GetVoltageDividerWeight(double RJunction, double R0, double Rground, unsigned z_index);
            double GetParallelResistance(std::vector<double> R, int NRESISTORS, int resistorindex);
        	virtual bool SetVoltageDampingFactors();
            virtual void SayHello();



    };


} /* namespace locust */

#endif /* LMCVOLTAGEDIVIDER_HH_ */
