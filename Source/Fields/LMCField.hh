/*
 * LMCField.hh
 *
 *  Created on: Jun. 4, 2021
 *      Author: pslocum
 */

#ifndef LMCFIELD_HH_
#define LMCFIELD_HH_

#include "param.hh"
#include "logger.hh"

#include <vector>

namespace locust
{
 /*!
 @class Field
 @author P. Slocum
 @brief Base class to characterize Field selection
 @details
 Available configuration options:
 No input parameters
 */


    class Field
    {

        public:
            Field();
            virtual ~Field();

            virtual bool Configure( const scarab::param_node& ){};

            virtual std::vector<double> TE_E(int l, int m, int n, double r, double theta, double z) const {};
            virtual std::vector<double> TE_H(int l, int m, int n, double r, double theta, double z) const {};
            virtual std::vector<double> TM_E(int l, int m, int n, double r, double theta, double z) const {};
            virtual std::vector<double> TM_H(int l, int m, int n, double r, double theta, double z) const {};
            virtual double Integrate(int l, int m, int n, bool teMode, bool eField){};

    };


}; /* namespace locust */

#endif /* LMCFIELD_HH_ */
