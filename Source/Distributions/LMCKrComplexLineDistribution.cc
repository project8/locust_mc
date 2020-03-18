/*
 * LMCKrComplexLineDistribution.cc
 *
 *  Created on: Mar 15, 2020
 *      Author: nbuzinsky
 */

#include "LMCKrComplexLineDistribution.hh"

#include "LMCConst.hh"

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/classification.hpp>

#include <iostream>

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
        fGases{"H2","Kr"},
        fEmittedPeak("shake"),
        fShakeInterpolator({0},{0},1)
        //fKrInterpolant(std::vector<double>(1).data(),std::vector<double>(1).data(),1,0),
    {
        if(aParam.has("fwhm"))
            fFWHM = aParam.get_value< double >( "fwhm", fFWHM );
        if(aParam.has("line-position"))
            fLinePosition = aParam.get_value< double >( "line-position", fLinePosition );
        if(aParam.has("emitted-peak"))
            fEmittedPeak = aParam.get_value< std::string >( "emitted-peak", fEmittedPeak );

        //create random number generator instances
        const double kr_line_width = 2.83; // eV
        fXArray = linspace(-1000,1000,fNPointsSELA);

        fNormal = std::normal_distribution<double>(0, gaussian_FWHM_to_sigma(fFWHM));
        fLorentzian = std::cauchy_distribution<double>(fLinePosition, kr_line_width / 2. ); //pos, hwhm
        fUniform = std::uniform_real_distribution<double>(0,1);

        for(unsigned i=0; i < fGases.size(); ++i)
            fGeometric[i] = std::geometric_distribution<int>(1. - fScatterProbability[i] );

        fShakeSpectrum = shake_spectrum();
        create_cdf(fShakeInterpolator, to_vector(fShakeSpectrum), to_vector(fXArray));

        for(unsigned i=0; i < fGases.size(); ++i)
        {
            fEnergyLossSpectra[i] = energy_loss_spectra(fGases[i]);
            //fix me (use correct (custom grid))
            create_cdf(fEnergyLossInterpolator[i], to_vector(fEnergyLossSpectra[i]), to_vector(fXArray));
        }


        for(unsigned i=0; i < fGases.size(); ++i)
            fGasIndex.insert( std::pair<std::string, unsigned>(fGases[i], i) );


    }


    void KrComplexLineDistribution::read_shake_data()
    {
        std::string filename = "KrShakeParameters214.txt";
        std::vector<std::vector<double> > read_data = read_file(filename, "," );

        //assign variables from these parameters
        //fGammaWidth

    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // The shakeup/shakeoff spectrum for the 17.8 keV line of Kr83m based on Hamish Vadantha paper
    // Overleaf link for the paper https://www.overleaf.com/project/5d56e015162d244b1c283c57
    // Read parameters for shake up shake off spectrum from an excel spread sheet.
    
    // nprime in Eq. (6)
    std::valarray<double> KrComplexLineDistribution::nprime(const double &E_b, const std::valarray<double> &W)
    {
        return sqrt(E_b/abs(W));
    }

    // Eq.(9) in Hamish Vedantha shake spectrum paper
    double KrComplexLineDistribution::C1s(const double &E_b)
    {
        return 98.2/E_b;
    }

    // Eq. (6)
    std::valarray<double> KrComplexLineDistribution::P_1s_nprime(const double &E_b, const std::valarray<double> &W)
    {
        std::valarray<double> n_prime = nprime(E_b, W);
        std::valarray<double> P_1s_nprime = C1s(E_b) / (1. - exp(-2. * LMCConst::Pi() * n_prime ));
        P_1s_nprime *= pow(n_prime,8);
        P_1s_nprime *= pow( pow(n_prime, 2) + 1., -4);
        P_1s_nprime *= exp(-4. * n_prime * atan(1. / n_prime));
        return P_1s_nprime;
    }

    // Eq. (5) in Hamish Vedantha shake spectrum paper
    // shake up spectrum for the ith state
    std::valarray<double> KrComplexLineDistribution::I(const unsigned &i, const std::valarray<double> &E)
    {
        double numerator = fAIntensity[i]*fGammaWidth[i];
        std::valarray<double> denominator = 2 * LMCConst::Pi() *( pow(fGammaWidth[i], 2.) /4. + pow(fECore - fBBinding[i] - E,2. ));
        return numerator/denominator;
    }

    // Eq. (11) in Hamish Vedantha shake spectrum paper
    // shake off spectrum for the ith state
    std::valarray<double> KrComplexLineDistribution::spectrum1(const unsigned &i, const std::valarray<double> &E)
    {
        std::valarray<double> spectrum;
        if (fEbScale[i] == 0 )
        {
            spectrum = I(i, E);
        }
        else
        {
            const double epsilon_shake_spectrum = 1e-4;
            std::valarray<double> factor = atan(2. * (fECore - E - fBBinding[i] + epsilon_shake_spectrum)/fGammaWidth[i]) / LMCConst::Pi() + 0.5;
            spectrum = fAIntensity[i] * factor * P_1s_nprime(fEbScale[i], fECore-E-fBBinding[i] + epsilon_shake_spectrum);
        }
        return spectrum;
    }

    // adding up shake spectrum for all the states
    std::valarray<double> KrComplexLineDistribution::full_shake_spectrum(const std::valarray<double> &E, const unsigned &start_number_of_i, const unsigned &end_number_of_i)
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

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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

    double KrComplexLineDistribution::EnergyLossSpectrum(double eLoss, double oscillator_strength)
    {
        const double T = fLinePosition; //approximate proxy instead of event by event
        return (eLoss > 0 ) ?  ((LMCConst::E_Rydberg() / eLoss) * oscillator_strength * log(4. * T * eLoss / pow(LMCConst::E_Rydberg(), 3.) )) : 0; // Produces energy loss spectrum
    }

    ////////////////////////////////////////////////////////////////////////////////////


    std::vector<std::vector<double> > KrComplexLineDistribution::read_file(std::string filename, std::string delimiter)
    {
        std::ifstream file(filename);
        std::vector<std::vector<double> > dataList;
        std::string line = "";
        // Iterate through each line and split the content using delimeter
        while (getline(file, line))
        {
            if(line.empty() || line.at(0) == '#' ) continue;
            std::vector<std::string> string_vec;
            boost::algorithm::split(string_vec, line, boost::is_any_of(delimiter));

            //change types
            std::vector<double> double_vec(string_vec.size());
            transform(string_vec.begin(), string_vec.end(), double_vec.begin(),
                    [](std::string const& val) {return stod(val);});
            dataList.push_back(double_vec);
        }

        file.close();
        return dataList;
    }

    std::valarray<double> KrComplexLineDistribution::energy_loss_spectra(const std::string &gas_species)
    {
        std::string filename = gas_species + "OscillatorStrength.txt";
        std::vector<std::vector<double> > read_data  = read_file(filename, "\t");
        std::sort(read_data.begin(), read_data.end(), [](const std::vector<double> & a, const std::vector<double> & b) -> bool { return a[0] < b[0]; });

    }


    double KrComplexLineDistribution::generate_from_cdf(double u, boost::math::barycentric_rational<double> &aCDF )
    {
        return aCDF(u);
    }
    

    // Reads in probability density, fills boost interpolator with corresponding inverse CDF for subsequent inversion sampling
    void KrComplexLineDistribution::create_cdf(boost::math::barycentric_rational<double> &interpolant, std::vector<double> f, std::vector<double> x)
    {
        std::vector<double> cdf = trapezoidal_rule(f,x);
        interpolant = boost::math::barycentric_rational<double>(cdf.data(), x.data(), x.size(), 0); //linear interpolation: don't change, doesnt converge
    }

    std::vector<double> KrComplexLineDistribution::trapezoidal_rule(std::vector<double> f, std::vector<double> x)
    {
        std::vector<double> cdf(x.size());

        for(int i=1;i<cdf.size();++i)
            cdf[i] = cdf[i-1] + (f[i-1] + f[i]) * (x[i] - x[i-1]) / 2.;

        double norm = cdf.back();

        std::transform(cdf.begin(), cdf.end(), cdf.begin(), [norm](double& c){return c/norm;});
        return cdf;
    }
    

    ////////////////////////////////////////////////////////////////////////////////////


    double KrComplexLineDistribution::generate_shake()
    {
        double u = fUniform(fRNEngine); //fix (why???)
        return generate_from_cdf(u,fShakeInterpolator);
    }

    std::string KrComplexLineDistribution::generate_gas_species()
    {
        double u = fUniform(fRNEngine);
        bool gas_choice = u < (fAmplitude[0] / (fAmplitude[0] + fAmplitude[1]));
        return fGases[int(gas_choice)];
    }

    double KrComplexLineDistribution::generate_energy_loss(std::string gas_species)
    {
        double u = fUniform(fRNEngine);
        unsigned gas_index = fGasIndex[gas_species];
        
        return generate_from_cdf(u, fEnergyLossInterpolator[gas_index]);
    }

    int KrComplexLineDistribution::generate_nscatters(std::string &gas_species)
    {
        return fGeometric[fGasIndex[gas_species]](fRNEngine);
    }

    double KrComplexLineDistribution::Generate()
    {
        //initialize energy from lineshape (lorentzian/ shake)
        double generated_energy;

        if(fEmittedPeak == "lorentzian")
            generated_energy = fLorentzian(fRNEngine);
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
        generated_energy += fNormal(fRNEngine);

        return generated_energy;

    }

    std::valarray<double> KrComplexLineDistribution::linspace(double a, double b, unsigned N)
    {
        double h = (b - a) / static_cast<double>(N-1);
        std::valarray<double> v(N);
        double val = a;
        for (unsigned i = 0; i < N; ++i)
        { 
            v[i] = val;
            val+=h;
        }
        return v;
    }

    // Converts a gaussian's FWHM to sigma
    double KrComplexLineDistribution::gaussian_FWHM_to_sigma(const double &fwhm)
    {
        return fwhm / (2. * sqrt(2. * log(2.)));
    }

    std::vector<double> KrComplexLineDistribution::to_vector(const std::valarray<double> a)
    {
        std::vector<double> v;
        v.assign(std::begin(a), std::end(a));
        return v;
    }


    ///////////////////////////////////////////////////////////////////////////////////////////

} /* namespace locust */
