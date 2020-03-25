/*
 * LMCRuddDistribution.hh
 *
 *  Created on: Mar 25, 2020
 *      Author: nbuzinsky
 */

#ifndef LMCRUDDDISTRIBUTION_HH_
#define LMCRUDDDISTRIBUTION_HH_

#include "LMCBaseDistribution.hh"

namespace locust
{
 /*!
 @class RuddDistribution
 @author N. Buzinsky
 @brief Angular distribution described in Eqn 14 of "Differential and total cross sections for ionization of helium and hydrogen by electrons", M.E. Rudd, Phys. Rev. A 44, 1644 – Published 1 August 1991. Useful for modeling angular distribution of electron-atom scattering
 @details
 Available configuration options:
 No input parameters
 */
    //template<class T>
    class RuddDistribution : public BaseDistribution
    {

        public:
            RuddDistribution(const scarab::param_node &aParam);
            double Generate();

        private:
            std::uniform_real_distribution<double> fDistribution;
            double fAlpha;
};


} /* namespace locust */

#endif /* LMCRUDDDISTRBUTION_HH_ */
