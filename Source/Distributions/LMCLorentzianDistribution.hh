/*
 * LMCLorentzianDistribution.hh
 *
 *  Created on: Mar 24, 2020
 *      Author: nbuzinsky
 */

#ifndef LMCLORENTZIANDISTRIBUTION_HH_
#define LMCLORENTZIANDISTRIBUTION_HH_

#include "LMCBaseDistribution.hh"

namespace locust
{
 /*!
 @class LorentzianDistribution
 @author N. Buzinsky
 @brief Derived class for Locust lorentzian/ cauchy distribution with location parameter "mean", scale parameter "fwhm"
 @details
 Available configuration options:
 No input parameters
 */
    //template<class T>
    class LorentzianDistribution : public BaseDistribution
    {

        public:
            //virtual T Generate() = 0;
            LorentzianDistribution(const scarab::param_node &aParam);
            double Generate();

        private:
            std::cauchy_distribution<double> fDistribution;
            double fMean;
            double fFWHM;
};


} /* namespace locust */

#endif /* LMCLORENTZIANDISTRBUTION_HH_ */
