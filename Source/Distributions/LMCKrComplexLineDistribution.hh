/*
 * LMCKrComplexLineDistribution.hh
 *
 *  Created on: Mar 12, 2020
 *      Author: nbuzinsky
 */

#ifndef LMCKRCOMPLEXLINEDISTRIBUTION_HH_
#define LMCKRCOMPLEXLINEDISTRIBUTION_HH_

#include "LMCBaseDistribution.hh"
#include "path.hh"
#include "macros.hh"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include <map>
#include <valarray>
#include <vector>

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
    class KrComplexLineDistribution : public BaseDistribution
    {

        public:
            KrComplexLineDistribution(const scarab::param_node &aParam);
            ~KrComplexLineDistribution();
            double Generate();

        private:
            double fFWHM;
            double fLinePosition;
            std::vector<double> fAmplitude;
            std::vector<double> fScatterProbability;
            std::string fEmittedPeak;

            unsigned fNPointsSELA;
            std::vector<std::string> fGases;
            std::map<std::string, unsigned> fGasIndex;

            std::string fDataDir;

            std::valarray<double> fXArray;
            std::valarray<double> fShakeSpectrum;
            std::vector<std::valarray<double> > fEnergyLossSpectra;

            std::vector<gsl_spline*> fEnergyLossInterpolator;
            std::vector<gsl_interp_accel*> fEnergyLossAccelerator;
            gsl_spline *fShakeInterpolator;
            gsl_interp_accel *fShakeAccelerator;
            
            std::vector<double> fAIntensity;
            std::vector<double> fBBinding;
            std::vector<double> fGammaWidth;
            std::vector<double> fEbScale;
            double fECore;


            std::uniform_real_distribution<double> fUniform;
            std::normal_distribution<double> fNormal;
            std::cauchy_distribution<double> fLorentzian;

            std::vector<std::geometric_distribution<int>> fGeometric;

            void read_shake_data();
            std::valarray<double> nprime(const double &E_b, const std::valarray<double> &W);
            double C1s(const double &E_b);
            std::valarray<double> P_1s_nprime(const double &E_b, const std::valarray<double> &W);
            std::valarray<double> I(const unsigned &i, const std::valarray<double> &E);
            std::valarray<double> spectrum1(const unsigned &i, const std::valarray<double> &E);
            std::valarray<double> full_shake_spectrum(const std::valarray<double> &E, const unsigned &start_number_of_i, const unsigned &end_number_of_i);
            std::valarray<double> shake_spectrum();

            double aseev_func_tail(const double &energy_loss, const std::string &gas_type);
            std::vector<double> EnergyLossSpectrum(std::vector<std::vector<double>> aData);
            std::vector<std::vector<double>> energy_loss_spectra(const std::string &gas_species);

            std::vector<std::vector<double> > read_file(std::string filename, std::string delimiter);
            double generate_from_cdf(double u, gsl_spline*& aCDFSpline, gsl_interp_accel*& aAccelerator );
            void create_cdf(gsl_spline*& interpolant, std::vector<double> f, std::vector<double> x);
            std::vector<double> trapezoidal_rule(std::vector<double> f, std::vector<double> x);

            double generate_shake();
            std::string generate_gas_species();
            double generate_energy_loss(std::string gas_species);
            int generate_nscatters(std::string &gas_species);

            std::vector<std::vector<double>> transpose_vector(const std::vector<std::vector<double>> aVector);
            std::valarray<double> linspace(double a, double b, unsigned N);
            double gaussian_FWHM_to_sigma(const double &fwhm);
            std::vector<double> to_vector(const std::valarray<double> a);
            void extrapolate_oscillator_strength(std::vector<std::vector<double>> &aOscillatorStrength, const std::string &gas_species);

};

} /* namespace locust */

#endif /* LMCKRCOMPLEXLINESHAPEDISTRBUTION_HH_ */
