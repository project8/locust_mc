/*
 * LMCUniformCircleDistribution.hh
 *
 *  Created on: Aug 12, 2020
 *      Author: nbuzinsky
 */

#ifndef LMCUNIFORMCIRCLEDISTRIBUTION_HH_
#define LMCUNIFORMCIRCLEDISTRIBUTION_HH_

#include "LMCBaseDistribution.hh"

namespace locust
{
 /*!
 @class UniformCircleDistribution
 @author N. Buzinsky
 @brief Derived class for Locust uniform circle distribution called between [0, radius]. Type must be double
 @details
 Available configuration options:
 No input parameters
 */
    //template<class T>
    class UniformCircleDistribution : public BaseDistribution
    {

        public:
            UniformCircleDistribution(const scarab::param_node &aParam);
            UniformCircleDistribution(const double &aRadius);
            double Generate();

        private:
            std::uniform_real_distribution<double> fDistribution;
            double fRadius;
};


} /* namespace locust */

#endif /* LMCUNIFORMCIRCLEDISTRBUTION_HH_ */
