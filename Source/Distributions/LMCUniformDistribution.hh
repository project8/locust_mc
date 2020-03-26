/*
 * LMCUniformDistribution.hh
 *
 *  Created on: Mar 10, 2020
 *      Author: nbuzinsky
 */

#ifndef LMCUNIFORMDISTRIBUTION_HH_
#define LMCUNIFORMDISTRIBUTION_HH_

#include "LMCBaseDistribution.hh"

namespace locust
{
 /*!
 @class UniformDistribution
 @author N. Buzinsky
 @brief Derived class for Locust uniform distribution called between [min_value, max_value]. Type can be double or int.
 @details
 Available configuration options:
 No input parameters
 */
    //template<class T>
    class UniformDistribution : public BaseDistribution
    {

        public:
            UniformDistribution(const scarab::param_node &aParam);
            UniformDistribution(const double &aMinValue, const double &aMaxValue);
            double Generate();

        private:
            std::uniform_real_distribution<double> fDistribution;
            double fMinValue;
            double fMaxValue;
};


} /* namespace locust */

#endif /* LMCUNIFORMDISTRBUTION_HH_ */
