/*
 * LMCDiracDistribution.hh
 *
 *  Created on: Mar 10, 2020
 *      Author: nbuzinsky
 */

#ifndef LMCDIRACDISTRIBUTION_HH_
#define LMCDIRACDISTRIBUTION_HH_

#include "LMCBaseDistribution.hh"

namespace locust
{
 /*!
 @class DiracDistribution
 @author N. Buzinsky
 @brief Derived class for Locust fixed distribution which trivially returns the fixed configured value from Generate()
 @details
 Available configuration options:
 No input parameters
 */
    //template<class T>
    class DiracDistribution : public BaseDistribution
    {

        public:
            DiracDistribution(const double &aValue);
            DiracDistribution(const scarab::param_node &aParam);
            double Generate();

        private:
            double fValue;
};


} /* namespace locust */

#endif /* LMCDIRACDISTRBUTION_HH_ */
