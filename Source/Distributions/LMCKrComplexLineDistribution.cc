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
    LOGGER( lmclog, "LMCKrComplexLineDistribution" );

    KrComplexLineDistribution::KrComplexLineDistribution(const scarab::param_node &aParam) :
        fECore( 0. ),
        fFWHM( 5. ),
        fLinePosition( 17826. ),
        fAmplitude{1,0},
        fScatterProbability(0.),
        fNPointsSELA(10000),
        fGases{"H2","Kr"},
        fEmittedPeak("shake")
    {
        if(aParam.has("fwhm"))
            fFWHM = aParam.get_value< double >( "fwhm", fFWHM );
        if(aParam.has("line-position"))
            fLinePosition = aParam.get_value< double >( "line-position", fLinePosition );
        if(aParam.has("emitted-peak"))
            fEmittedPeak = aParam.get_value< std::string >( "emitted-peak", fEmittedPeak );
        if(aParam.has("scatter-probability"))
            fScatterProbability = aParam.get_value< double >( "scatter-probability", fScatterProbability );

        for(unsigned i=0; i < fGases.size(); ++i)
            fGasIndex.insert( std::pair<std::string, unsigned>(fGases[i], i) );

        if(aParam.has("amplitude"))
        {
            scarab::param_node aAmplitudeNode = aParam["amplitude"].as_node();
            for(unsigned i=0;i<fGases.size();++i)
                if(aAmplitudeNode.has(fGases[i]))
                    fAmplitude[fGasIndex[fGases[i]]] = aAmplitudeNode.get_value<double >( fGases[i], fAmplitude[i]) ;
        }

        //create random number generator instances
        const double kr_line_width = 2.83; // eV
        fXArray = linspace(-1000,1000,fNPointsSELA);

        fNormal = std::normal_distribution<double>(0, gaussian_FWHM_to_sigma(fFWHM));
        fLorentzian = std::cauchy_distribution<double>(0, kr_line_width / 2. ); //pos, hwhm
        fUniform = std::uniform_real_distribution<double>(0,1);

        fGeometric = std::geometric_distribution<int>(1. - fScatterProbability);

        fShakeAccelerator = gsl_interp_accel_alloc();

        fDataDir =  TOSTRING(PB_DATA_INSTALL_DIR);

        //initialize shakeon/ shakeoff data
        read_shake_data();
        fShakeSpectrum = shake_spectrum();
        create_cdf(fShakeInterpolator, to_vector(fShakeSpectrum), to_vector(fXArray));

        //initialize energy loss scattering data
        fEnergyLossInterpolator = std::vector<gsl_spline*>(fGases.size());
        fEnergyLossAccelerator = std::vector<gsl_interp_accel*>(fGases.size());

        for(unsigned i=0; i < fGases.size(); ++i)
        {
            std::vector<std::vector<double> > scattering_data = energy_loss_spectra(fGases[i]);
            create_cdf(fEnergyLossInterpolator[i], scattering_data[1], scattering_data[0]);
        }

        std::string amp_log_string, scatter_log_string;
        for(unsigned i=0;i<fGases.size();++i)
        {
            amp_log_string += std::to_string(fAmplitude[i]) + "  ";
        }

        LDEBUG( lmclog, "Created kr-complex-line distribution. fwhm: " <<fFWHM<<" line-position: "<<fLinePosition<<" emitted-peak: "<<fEmittedPeak<<" scatter-prob: "<<fScatterProbability);
        LDEBUG( lmclog, "amplitudes: {" <<amp_log_string<<"} scatter-probabilities: {"<<scatter_log_string<<"}");
    }

    KrComplexLineDistribution::~KrComplexLineDistribution()
    {
        for(unsigned i=0;i<fGases.size(); ++i)
        {
            gsl_spline_free(fEnergyLossInterpolator[i]);
            gsl_interp_accel_free(fEnergyLossAccelerator[i]);
        }

        gsl_spline_free(fShakeInterpolator);
        gsl_interp_accel_free(fShakeAccelerator);
    }

    void KrComplexLineDistribution::read_shake_data()
    {
        std::string filename = fDataDir+"/KrShakeParameters214.txt";
        std::vector<std::vector<double> > read_data = read_file(filename, "," );
        read_data = transpose_vector(read_data);

        //assign variables from shakeoff data file
        fAIntensity = read_data[1];
        fBBinding = read_data[2];
        fEbScale = read_data[3];
        fGammaWidth = read_data[4];
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
        std::valarray<double> ones(1., P_1s_nprime.size());

        P_1s_nprime *= pow(n_prime, std::valarray<double>(8. * ones));
        P_1s_nprime *= pow( pow(n_prime, std::valarray<double>(2. * ones)) + 1., std::valarray<double>(-4. * ones));
        P_1s_nprime *= exp(-4. * n_prime * atan(1. / n_prime));
        return P_1s_nprime;
    }

    // Eq. (5) in Hamish Vedantha shake spectrum paper
    // shake up spectrum for the ith state
    std::valarray<double> KrComplexLineDistribution::I(const unsigned &i, const std::valarray<double> &E)
    {
        double numerator = fAIntensity[i]*fGammaWidth[i];
        std::valarray<double> ones(1., E.size());
        std::valarray<double> denominator = 2 * LMCConst::Pi() *( pow(fGammaWidth[i], std::valarray<double>(2. * ones)) /4. + pow(fECore - fBBinding[i] - E, std::valarray<double>(2. * ones) ));
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
        std::valarray<double> full_spectrum(double(0), E.size());
        for( unsigned i = start_number_of_i; i < end_number_of_i; ++i)
            full_spectrum += spectrum1(i, E);
        return full_spectrum;
    }

    // shake spectrum by adding up shake spectrum for all the states up to i=24
    std::valarray<double> KrComplexLineDistribution::shake_spectrum()
    {
        return full_shake_spectrum(fXArray, 0, 24);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //energy loss code

    // A sub function for the scatter function. Found in 
    // "Energy loss of 18 keV electrons in gaseous T and quench condensed D films" 
    // by V.N. Aseev et al. 2000
    double KrComplexLineDistribution::aseev_func_tail(const double &energy_loss, const std::string & gas_type)
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
        return A2*pow(omeg2, 2.)/(pow(omeg2, 2.) + 4. * pow(energy_loss-eps2, 2.) );
    }

    std::vector<double> KrComplexLineDistribution::EnergyLossSpectrum(std::vector<std::vector<double>> aData)
    {
        const double T = fLinePosition; //approximate proxy instead of event by event
        std::vector<double> aSpectrum;
        for(unsigned i=0; i<aData[0].size();++i)
        {
            double eLoss = aData[0][i];
            if(eLoss>0)
            {
                double oscillator_strength = aData[1][i];
                aSpectrum.push_back((LMCConst::E_Rydberg() / eLoss) * fabs(oscillator_strength) * log(4. * T * eLoss / pow(LMCConst::E_Rydberg(), 3.) )); // Produces energy loss spectrum
            }
            else
            {
                aSpectrum.push_back(0);
            }
        }
        return aSpectrum;
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

    std::vector<std::vector<double>> KrComplexLineDistribution::energy_loss_spectra(const std::string &gas_species)
    {
        std::string filename = fDataDir + "/" + gas_species + "OscillatorStrength.txt";
        std::vector<std::vector<double> > read_data  = read_file(filename, "\t");

        std::sort(read_data.begin(), read_data.end(), [](const std::vector<double> & a, const std::vector<double> & b) -> bool { return a[0] < b[0]; });

        //remove duplicates
        auto equal_lambda = [](const std::vector<double> &a, const std::vector<double> & b) { return a[0] == b[0];};
        read_data.erase( std::unique( read_data.begin(), read_data.end(), equal_lambda ), read_data.end() );

        read_data = transpose_vector(read_data);

        extrapolate_oscillator_strength(read_data, gas_species );

        read_data[1] = EnergyLossSpectrum(read_data);
        return read_data;

    }


    void KrComplexLineDistribution::extrapolate_oscillator_strength(std::vector<std::vector<double>> &aOscillatorStrength, const std::string &gas_species)
    {
        double max_energy = aOscillatorStrength[0].back();
        const double dE = 1;
        while( max_energy < 1000)
        {
            max_energy  += dE;
            aOscillatorStrength[0].push_back(max_energy);
            aOscillatorStrength[1].push_back(aseev_func_tail(max_energy, gas_species));
        }
        return;

    }


    double KrComplexLineDistribution::generate_from_cdf(double u, gsl_spline*& aCDFSpline, gsl_interp_accel*& aAccelerator )
    {
        return gsl_spline_eval(aCDFSpline, u, aAccelerator);
    }
    

    // Reads in probability density, fills boost interpolator with corresponding inverse CDF for subsequent inversion sampling
    void KrComplexLineDistribution::create_cdf(gsl_spline*& interpolant, std::vector<double> f, std::vector<double> x)
    {
        //ensure pdf is strictly increasing
        for(auto it=f.begin();it<f.end();++it)
        {
            if(*it <= 0)
            {
                unsigned ind = it - f.begin();
                f.erase(f.begin() + ind);
                x.erase(x.begin() + ind);
            }
        }

        std::vector<double> cdf = trapezoidal_rule(f,x);
        interpolant = gsl_spline_alloc(gsl_interp_cspline, cdf.size());
        gsl_spline_init(interpolant, cdf.data(), x.data(), x.size());
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
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    double KrComplexLineDistribution::generate_shake()
    {
        double u = fUniform(*fRNEngine);
        return generate_from_cdf(u,fShakeInterpolator, fShakeAccelerator);
    }

    std::string KrComplexLineDistribution::generate_gas_species()
    {
        double u = fUniform(*fRNEngine);
        bool gas_choice = u > (fAmplitude[0] / (fAmplitude[0] + fAmplitude[1]));
        return fGases[int(gas_choice)];
    }

    double KrComplexLineDistribution::generate_energy_loss(std::string gas_species)
    {
        double u = fUniform(*fRNEngine);
        unsigned gas_index = fGasIndex[gas_species];
        
        return generate_from_cdf(u, fEnergyLossInterpolator[gas_index], fEnergyLossAccelerator[gas_index]);
    }

    int KrComplexLineDistribution::generate_nscatters()
    {
        return fGeometric(*fRNEngine);
    }

    double KrComplexLineDistribution::Generate()
    {
        //initialize energy from lineshape (lorentzian/ shake)
        double generated_energy = fLinePosition;

        if(fEmittedPeak == "lorentzian")
            generated_energy += fLorentzian(*fRNEngine);
        else if(fEmittedPeak == "shake")
            generated_energy += generate_shake();
        else if(fEmittedPeak == "dirac" || fEmittedPeak == "fixed") //just to be explicit
            generated_energy += 0;

        //choose gas species
        std::string gas_species = generate_gas_species();
        
        //generate number of scatters
        int nScatters = generate_nscatters();

        for(int i=0; i < nScatters; ++i)
        {
            gas_species = generate_gas_species();
            generated_energy -= generate_energy_loss(gas_species);
        }

        //include gaussian smearing from finite detector resolution
        generated_energy += fNormal(*fRNEngine);

        return generated_energy;

    }

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::vector<std::vector<double>> KrComplexLineDistribution::transpose_vector(const std::vector<std::vector<double>> aVector)
    {
        std::vector<std::vector<double>> aArrays;
        for(unsigned i=0; i<aVector[0].size(); ++i)
        {
            std::vector<double> v;
            for(unsigned j=0; j<aVector.size(); ++j)
                v.push_back(aVector[j][i]);

            aArrays.push_back(v);
        }
        return aArrays;
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
