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

            constexpr double nprime(const double &E_b, const double &W);
            constexpr double C1s(const double &E_b);
            double P_1s_nprime(const double &E_b, const double &W);
            double I(const unsigned &i, const double &E);
            double spectrum1(const unsigned &i, const double &E);
            std::valarray<double> full_shake_spectrum(const double &E, const unsigned &start_number_of_i, const unsigned &end_number_of_i);
            std::valarray<double> shake_spectrum();
            std::valarray<double> aseev_func_tail(std::valarray<double> energy_loss_array, std::string gas_type);
            double EnergyLossSpectrum(double eLoss, double oscillator_strength);
            std::vector<std::vector<std::double> > read_file(std::string filename, std::string delimiter);
            std::valarray<double> energy_loss_spectra(const std::string &gas_species);
            double generate_from_cdf(double u, boost::math::barycentric_rational<double> &aCDF );
            void create_cdf(boost::math::barycentric_rational<double> &interpolant, std::vector<double> f, std::vector<double> x);

            template <typename T>
            T<double> trapezoidal_rule(T<double> f, T<double> x);

            double generate_shake();
            std::string generate_gas_species();
            double generate_energy_loss(std::string gas_species);
            int generate_nscatters(std::string &gas_species);

            template <typename T>
            std::vector<T> linspace(T a, T b, unsigned N);

            constexpr double gaussian_FWHM_to_sigma(const double &fwhm);



};


} /* namespace locust */

#endif /* LMCUNIFORMDISTRBUTION_HH_ */
