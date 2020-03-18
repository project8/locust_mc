/*
 * LMCKrComplexLineDistribution.cc
 *
 *  Created on: Mar 15, 2020
 *      Author: nbuzinsky
 */

#include "LMCKrComplexLineDistribution.hh"

#include <iostream>
#include <valarray>
#include <vector>

namespace locust
{

    //ReadFile((dataDir / "KrOscillatorStrength.txt").string(), krData);

    //pass gas parameters as arrays?
    KrComplexLineDistribution::KrComplexLineDistribution(const scarab::param_node &aParam) :
        fFWHM( 5. ),
        fLinePosition( 17826. ),
        fAmplitude(),
        fScatterProbability(),
        fNPointsSELA(10000),
        fGases("H2","Kr"),
        fEmittedPeak("shake")
    {
        if(aParam.has("fwhm"))
            fFWHM = aParam.get_value< double >( "fwhm", fFWHM );
        if(aParam.has("line-position"))
            fLinePosition = aParam.get_value< double >( "line-position", fLinePosition );
        if(aParam.has("emitted-peak"))
            fEmittedPeak = aParam.get_value< std::string >( "emitted-peak", fEmittedPeak );

        //create rng instances
        const double kr_line_width = 2.83; // eV
        fXArray = linspace(-1000,1000,fNPointsSELA);

        fNormal = std::normal_distribution<double>(0, gaussian_FWHM_to_sigma(fFWHM));
        fLorentzian = std::cauchy_distribution<double>(fLinePosition, kr_line_width / 2. ); //pos, hwhm
        fGeometric = std::geometric_distribution<int>(1. - fScatterProbability );
        fUniform = std::uniform_real_distribution<double>(0,1);

        fShakeSpectrum = shake_spectrum();


    }

    // The shakeup/shakeoff spectrum for the 17.8 keV line of Kr83m based on Hamish Vadantha paper
    // Overleaf link for the paper https://www.overleaf.com/project/5d56e015162d244b1c283c57
    // Read parameters for shake up shake off spectrum from an excel spread sheet.
    
    // nprime in Eq. (6)
    constexpr double KrComplexLineDistribution::nprime(const double &E_b, const double &W)
    {
        return sqrt(E_b/fabs(W));
    }

    // Eq.(9) in Hamish Vedantha shake spectrum paper
    constexpr double KrComplexLineDistribution::C1s(const double &E_b)
    {
        return 98.2/E_b;
    }

    // Eq. (6)
    double KrComplexLineDistribution::P_1s_nprime(const double &E_b, const double &W)
    {
        double n_prime = nprime(E_b, W);
        double P_1s_nprime = 1;
        P_1s_nprime *= C1s(E_b);
        P_1s_nprime /= (1. - exp(-2. * LMCConst::Pi() * n_prime ));
        P_1s_nprime *= pow(n_prime,8);
        P_1s_nprime *= pow( pow(n_prime, 2) + 1., -4);
        P_1s_nprime *= exp(-4. * n_prime * atan(1. / n_prime));
        return P_1s_nprime;
    }

    // Eq. (5) in Hamish Vedantha shake spectrum paper
    // shake up spectrum for the ith state
    double KrComplexLineDistribution::I(const unsigned &i, const double &E)
    {
        double numerator = A_Intensity[i]*Gamma_Width[i];
        double denominator = 2 * LMCConst::Pi() *( pow(Gamma_Width[i], 2.) /4. + pow(Ecore - B_Binding[i] - E,2. ));
        return numerator/denominator;
    }

    // Eq. (11) in Hamish Vedantha shake spectrum paper
    // shake off spectrum for the ith state
    double KrComplexLineDistribution::spectrum1(const unsigned &i, const double &E)
    {
        double spectrum;
        if (E_b_Scale[i] == 0 )
        {
            spectrum = I(i, E);
        }
        else
        {
            double factor = atan(2. * (Ecore - E - B_Binding[i] + epsilon_shake_spectrum)/Gamma_Width[i]) / LMCConst::Pi() + 0.5;
            spectrum = A_Intensity[i] * factor * P_1s_nprime(E_b_Scale[i], Ecore-E-B_Binding[i] + epsilon_shake_spectrum);
        }
        return spectrum;
    }

    // adding up shake spectrum for all the states
    std::valarray<double> KrComplexLineDistribution::full_shake_spectrum(const double &E, const unsigned &start_number_of_i, const unsigned &end_number_of_i)
    {
        std::valarray<double> full_spectrum;
        for( unsigned i = start_number_of_i; i < end_number_of_i; ++i)
            full_spectrum += spectrum1(i, E);
        return full_spectrum;
    }

    // shake spectrum by adding up shake spectrum for all the states up to i=24
    std::valarray<double> KrComplexLineDistribution::shake_spectrum()
    {
        //x_array = flip_array(x_array);
        return full_shake_spectrum(fXArray, 0, 24);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////
    //energy loss code

    // A sub function for the scatter function. Found in 
    // "Energy loss of 18 keV electrons in gaseous T and quench condensed D films" 
    // by V.N. Aseev et al. 2000
    std::valarray<double> KrComplexLineDistribution::aseev_func_tail(std::valarray<double> energy_loss_array, std::string gas_type)
    {
        double A2, omeg2, eps2;
        if(gas_type == "H2")
        {
            A2  = 0.195; omeg2 = 14.13; eps2 = 10.60;
        }
        else if(gas_type=="Kr")
        {
            A2 = 0.4019; omeg2 = 22.31; eps2 = 16.725;
        }
        return A2*pow(omeg2, 2.)/(pow(omeg2, 2.) + 4. * pow(energy_loss_array-eps2, 2.) );
    }

    //should keep eLoss > 0 fix
    double KrComplexLineDistribution::EnergyLossSpectrum(double eLoss, double oscillator_strength)
    {
        double T = rel_energy(fStartFrequencyMax, fBField);
        return (eLoss > 0 ) ?  ((LMCConst::E_Rydberg() / eLoss) * oscillator_strength * log(4. * T * eLoss / pow(LMCConst::E_Rydberg(), 3.) )) : 0; // Produces energy loss spectrum
    }

    ////////////////////////////////////////////////////////////////////////////////////

        //std::sort(readData.begin(), readData.end(), [](const std::pair<double,double> & a, const std::pair<double, double> & b) -> bool { return a.first < b.first; });

    std::vector<std::vector<std::string> > KrComplexLineDistribution::read_file(std::string filename, std::string delimiter)
    {
        std::ifstream file(filename);
        std::vector<std::vector<std::string> > dataList;
        std::string line = "";
        // Iterate through each line and split the content using delimeter
        while (getline(file, line))
        {
            if(line.empty() || line[0] == std::string("#")) continue;
            std::vector<std::string> vec;
            boost::algorithm::split(vec, line, boost::is_any_of(delimeter));
            dataList.push_back(vec);
        }

        file.close();
        return dataList;
    }


    double KrComplexLineDistribution::generate_from_cdf(double u, boost::math::barycentric_rational<double> &aCDF )
    {
        return aCDF(u);
    }
    

    // Reads in probability density, fills boost interpolator with corresponding inverse CDF for subsequent inversion sampling
    void KrComplexLineDistribution::create_cdf(boost::math::barycentric_rational<double> &interpolant, std::vector<double> f, std::vector<double> x)
    {
        std::vector<double> cdf(x.size());
        //double normalization = trapezoidal_rule(f,x);


        interpolant = boost::math::barycentric_rational<double>(cdf.data(), energies.data(), energies.size(), 0); //linear interpolation: don't change, doesnt converge
    }

    template <typename T>
    T<double> KrComplexLineDistribution::trapezoidal_rule(T<double> f, T<double> x)
    {
        T<double> cdf(x.size());

        for(int i=1;i<cdf.size();++i)
            cdf[i] = cdf[i-1] + (energy_loss[i-1] + energy_loss[i]) * (x[i] - x[i-1]) / 2.;

        double norm = cdf.back();

        std::transform(cdf.begin(), cdf.end(), cdf.begin(), [norm](double& c){return c/norm;});
        return cdf;
    }
    

    ////////////////////////////////////////////////////////////////////////////////////


    double KrComplexLineDistribution::generate_shake()
    {
         //return fLorentizian(fGenerator);
    }

    std::string KrComplexLineDistribution::generate_gas_species()
    {
        double u = fUniform(fGenerator);
        bool gas_choice = u < (fAmplitude[0] / (fAmplitude[0] + fAmplitude[1]));
        return fGases[int(gas_choice)];
    }

    double generate_energy_loss(std::string gas_species)
    {
        double u = fUniform(fGenerator);
        unsigned gas_index = gas_index[gas_species];
        
        return generate_from_cdf(u, fEnergyLossInterpolator[gas_index]);
    }

    int KrComplexLineDistribution::generate_nscatters(std::string &gas_species)
    {
        return fGeometric[species_index[gas_species]](fGenerator);
    }

    double KrComplexLineDistribution::Generate()
    {
        //initialize energy from lineshape (lorentzian/ shake)
        double generated_energy;

        if(fEmittedPeak == "lorentzian")
            generated_energy = fLorentizian(fGenerator);
        if(fEmittedPeak == "shake")
            generated_energy = generate_shake();
        
        
        //choose gas species
        std::string gas_species = generate_gas_species();
        
        //generate number of scatters
        int nScatters = generate_nscatters(gas_species);

        //calculate n missed tracks on energy loss
        for(int i=0; i < nScatters; ++i)
            generated_energy -= generate_energy_loss(gas_species);
        
        //include gaussian smearing from finite detector resolution
        generated_energy += fNormal(fGenerator);

        return generated_energy;

    }

    template <typename T>
    std::vector<T> KrComplexLineDistribution::linspace(T a, T b, unsigned N)
    {
        T h = (b - a) / static_cast<T>(N-1);
        std::vector<T> v(N);
        for (unsigned i = 0, T val = a; i < N; ++i, val += h) v[i] = val;
        return v;
    }

    // Converts a gaussian's FWHM to sigma
    constexpr double KrComplexLineDistribution::gaussian_FWHM_to_sigma(const double &fwhm)
    {
        return fwhm / (2. * sqrt(2. * log(2.)));
    }


    ///////////////////////////////////////////////////////////////////////////////////////////

} /* namespace locust */
