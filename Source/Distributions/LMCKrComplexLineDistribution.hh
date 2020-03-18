/*
 * LMCKrComplexLineDistribution.hh
 *
 *  Created on: Mar 12, 2020
 *      Author: nbuzinsky
 */

#ifndef LMCKRCOMPLEXLINEDISTRIBUTION_HH_
#define LMCKRCOMPLEXLINEDISTRIBUTION_HH_

#include "LMCBaseDistribution.hh"

namespace locust
{
 /*!
 @class KrComplexLineDistribution
 @author N. Buzinsky
 @brief Derived class for custom P8 (Phase II) specific complex line shape using Kr
 @details
 Available configuration options:
 No input parameters
 */
    //template<class T>
    //////change fGenerator name!!!!!???
    class KrComplexLineDistribution : public BaseDistribution
    {

        public:
            KrComplexLineDistribution(const scarab::param_node &aParam);
            double Generate();

        private:
            double fFWHM;
            double fLinePosition;
            std::vector<double> fAmplitude;
            std::vector<double> fScatterProbability;
            unsigned fNPointsSELA;
            std::vector<std::string> fGases;

            std::vector<boost::math::barycentric_rational<double> > fEnergyLossInterpolator;
            boost::math::barycentric_rational<double> fShakeInterpolator;
            
            std::string fEmittedPeak;

            std::valarray<double> fXArray;
            std::valarray<double> fShakeSpectrum;
            std::vector<std::valarray<double> > fEnergyLossSpectra;


            std::uniform_real_distribution<double> fUniform;
            std::normal_distribution<double> fNormal;
            std::cauchy_distribution<double> fLorentzian;

            std::vector<std::normal_distribution<double>> fGeometric;
};


} /* namespace locust */

#endif /* LMCUNIFORMDISTRBUTION_HH_ */
