/*
 * LMCFakeTrackSignalGenerator.cc
 *
 *  Created on: Aug 8 2018
 *      Author: plslocum, buzinsky
 */

#include "LMCFakeTrackSignalGenerator.hh"
#include "LMCDigitizer.hh"
#include "logger.hh"
#include "LMCConst.hh"
#include "LMCThreeVector.hh"
#include "path.hh"
#include <algorithm>
#include <random>
#include <math.h>
#include <sstream>
#include <string>
#include <iostream>
#include <numeric>

using std::string;

namespace locust
{
    LOGGER( lmclog, "FakeTrackSignalGenerator" );

    MT_REGISTER_GENERATOR(FakeTrackSignalGenerator, "fake-track");

    FakeTrackSignalGenerator::FakeTrackSignalGenerator( const std::string& aName ) :
        Generator( aName ),
        fDoGenerateFunc( &FakeTrackSignalGenerator::DoGenerateTime ),
        fSignalPower( 0. ),
        fStartVPhase( 0. ),
        fStartTimeMin( 0. ),
        fStartTimeMax( 0. ),
        fStartPitchMin( 89.9 ),
        fStartPitchMax( 90. ),
        fPitchMin( 0. ),
        fLO_frequency( 0. ),
        fNTracksMean( 0. ),
        fBField(1.0),
        fRandomSeed(0),
        fNEvents(1),
        fAharmonicCorrection( false ),
        fAharmonicPowerCoupling( 1. ),
        fPitchCorrection( true ),
        fSlopeCorrection( false ),
        fRandomEngine(0),
        fHydrogenFraction(1),
        fTrapLength(0.1784),  //Phase II harmonic trap L0 (A. Ashtari Esfahani et al.- Phys. Rev. C 99, 055501 )
        fRootFilename("LocustEvent.root"),
        fSlope( 0. ),
        fPitch( 0. ),
        fTrackLength( 0. ),
        fStartTime( 0. ),
        fEndTime( 0. ),
        fStartFrequency( 0. ),
        fCoilCurrent( 0.3 ),
        fCurrentFrequency( 0. ),
        fUseEnergyDistribution(false),
        fUseFrequencyDistribution(false),
        fNTracks(0)
    {
        fRequiredSignalState = Signal::kTime;
    }

    FakeTrackSignalGenerator::~FakeTrackSignalGenerator()
    {
        for(unsigned i=0;i<fInterpolators.size(); ++i)
        {
            gsl_spline_free(fInterpolators[i]);
            gsl_interp_accel_free(fAccelerators[i]);
        }
    }

    bool FakeTrackSignalGenerator::Configure( const scarab::param_node& aParam )
    {
        if( aParam.has( "scattering-angle" ) )
        {
            fScatteringAngleDistribution = fDistributionInterface.get_dist(aParam["scattering-angle"].as_node());
        }
        else
        {
            LWARN( lmclog, "Using default distribution: Scattering Angle = 0 ");
            scarab::param_node default_setting;
            default_setting.add("name","dirac");
            fScatteringAngleDistribution = fDistributionInterface.get_dist(default_setting);
        }

        if( aParam.has( "start-frequency" ) )
        {
            fStartFrequencyDistribution = fDistributionInterface.get_dist(aParam["start-frequency"].as_node());
            fUseFrequencyDistribution = true;
        }
        if( aParam.has( "start-energy" ) )
        {
            fStartEnergyDistribution = fDistributionInterface.get_dist(aParam["start-energy"].as_node());
            fUseEnergyDistribution = true;
        }

        if( aParam.has( "slope" ) )
        {
            fSlopeDistribution = fDistributionInterface.get_dist(aParam["slope"].as_node());
        }
        else
        {
            LWARN( lmclog, "Using default distribution: Slope = 0 ");
            scarab::param_node default_setting;
            default_setting.add("name","dirac");
            fSlopeDistribution = fDistributionInterface.get_dist(default_setting);
        }

        if( aParam.has( "track-length" ) )
        {
            fTrackLengthDistribution = fDistributionInterface.get_dist(aParam["track-length"].as_node());
        }
        else
        {
            LWARN( lmclog, "Using default distribution: Track Length = 1e-4 ");
            scarab::param_node default_setting;
            default_setting.add("name","dirac");
            default_setting.add("value","1e-4");
            fTrackLengthDistribution = fDistributionInterface.get_dist(default_setting);
        }

        if( aParam.has( "radius" ) )
        {
            fRadiusDistribution = fDistributionInterface.get_dist(aParam["radius"].as_node());
        }
        else
        {
            LWARN( lmclog, "Using default distribution: radius = 0 ");
            scarab::param_node default_setting;
            default_setting.add("name","dirac");
            fRadiusDistribution = fDistributionInterface.get_dist(default_setting);
        }

        if( aParam.has( "z0" ) )
        {
            fz0Distribution = fDistributionInterface.get_dist(aParam["z0"].as_node());
        }
        else
        {
            LWARN( lmclog, "Using default distribution: z0 = 0 ");
            scarab::param_node default_setting;
            default_setting.add("name","dirac");
            fz0Distribution = fDistributionInterface.get_dist(default_setting);
        }

        if( aParam.has( "signal-power" ) )
            SetSignalPower( aParam.get_value< double >( "signal-power", fSignalPower ) );

        if( aParam.has( "start-vphase" ) )
            SetStartVPhase( aParam.get_value< double >( "start-vphase", fStartVPhase ) );

        if( aParam.has( "start-pitch-max" ) )
            SetStartPitchMax( aParam.get_value< double >( "start-pitch-max", fStartPitchMax ) );

        if( aParam.has( "start-pitch-min" ) )
            SetStartPitchMin( aParam.get_value< double >( "start-pitch-min", fStartPitchMin ) );

        if( aParam.has( "min-pitch" ) )
            SetPitchMin( aParam.get_value< double >( "min-pitch", fPitchMin ) );

        if( aParam.has( "start-time-max" ) )
            SetStartTimeMax( aParam.get_value< double >( "start-time-max", fStartTimeMax ) );

        if( aParam.has( "start-time-min" ) )
            SetStartTimeMin( aParam.get_value< double >( "start-time-min", fStartTimeMin ) );

        if( aParam.has( "lo-frequency" ) )
            SetFrequency( aParam.get_value< double >( "lo-frequency", fLO_frequency ) );

        if (aParam.has( "ntracks-mean") )
            SetNTracksMean( aParam.get_value< double >( "ntracks-mean",fNTracksMean) );

        if (aParam.has("magnetic-field") )
            SetBField(  aParam.get_value< double >("magnetic-field", fBField) );

        if (aParam.has( "random-seed") )
            SetRandomSeed(  aParam.get_value< int >( "random-seed",fRandomSeed) );

        if (aParam.has( "n-events") )
            SetNEvents(  aParam.get_value< int >( "n-events",fNEvents) );

        if (aParam.has( "aharmonic-correction") )
            SetAharmonicCorrection(  aParam.get_value< bool >( "aharmonic-correction", fAharmonicCorrection) );

        if (aParam.has( "pitch-correction") )
            SetPitchCorrection(  aParam.get_value< bool >( "pitch-correction", fPitchCorrection) );

        if (aParam.has( "slope-correction") )
            SetSlopeCorrection( aParam.get_value< bool >( "slope-correction",fSlopeCorrection) );

        if (aParam.has( "coil-current") )
            SetCoilCurrent(  aParam.get_value< double >( "coil-current", fCoilCurrent) );

        if (aParam.has( "trap-length") )
            SetTrapLength(  aParam.get_value< double >( "trap-length", fTrapLength) );

        if (aParam.has( "hydrogen-fraction") )
        {
            SetHydrogenFraction(  aParam.get_value< int >( "hydrogen-fraction",fHydrogenFraction) );

            if( fHydrogenFraction > 1 ||  fHydrogenFraction < 0)
                LERROR( lmclog, "hydrogen-fraction must be between 0 and 1!");
        }

        if( aParam.has( "root-filename" ) )
        {
            fRootFilename = aParam["root-filename"]().as_string();
        }

        if(fRandomSeed)
            fRandomEngine.seed(fRandomSeed);
        else
        {
            std::random_device rd;
            fRandomEngine.seed(rd());
        }

        fDistributionInterface.SetSeed(fRandomSeed);

        if(!fNTracksMean && !fPitchMin)
            LERROR( lmclog, "No condition set for NTracks per event! Set one of pitch-min or ntracks-mean");

        if( fUseFrequencyDistribution && fUseEnergyDistribution)
            LERROR( lmclog, "User specified both start frequency and start energy distribution! Please specify only one!");

        if( ! (fUseFrequencyDistribution || fUseEnergyDistribution))
        {
            LWARN( lmclog, "Using default distribution: Frequency = fLO + 50 MHz ");
            scarab::param_node default_setting;
            default_setting.add("name","dirac");
            default_setting.add("value",fLO_frequency + 50e6);
            fStartFrequencyDistribution = fDistributionInterface.get_dist(default_setting);
            fUseFrequencyDistribution = true;
        }

        if( aParam.has( "domain" ) )
        {
            string domain = aParam["domain"]().as_string();
            if( domain == "time" )
            {
                SetDomain( Signal::kTime );
                LDEBUG( lmclog, "Domain is equal to time.");
            }
            else if( domain == "freq" )
            {
                SetDomain( Signal::kFreq );
            }
            else
            {
                LERROR( lmclog, "Unable to use domain requested: <" << domain << ">" );
                return false;
            }
        }

        std::vector<std::pair<double, double> > h2Data, krData;
        scarab::path dataDir = aParam.get_value( "data-dir", ( TOSTRING(PB_DATA_INSTALL_DIR) ) );
        LDEBUG( lmclog, "Data directory: " << dataDir );
        ReadFile((dataDir / "H2OscillatorStrength.txt").string(), h2Data);
        ReadFile((dataDir / "KrOscillatorStrength.txt").string(), krData);
        ExtrapolateData(h2Data, std::array<double, 3>{0.195, 14.13, 10.60});
        ExtrapolateData(krData, std::array<double, 3>{0.4019, 22.31, 16.725});
        fInterpolators = std::vector<gsl_spline*>(2);
        fAccelerators = std::vector<gsl_interp_accel*>(2);
        SetInterpolator(fInterpolators[0],h2Data);
        SetInterpolator(fInterpolators[1],krData);

        return true;
    }

    void FakeTrackSignalGenerator::Accept( GeneratorVisitor* aVisitor ) const
    {
        aVisitor->Visit( this );
        return;
    }

    double FakeTrackSignalGenerator::GetSignalPower() const
    {
        return fSignalPower;
    }

    void FakeTrackSignalGenerator::SetSignalPower( double aPower )
    {
        fSignalPower = aPower;
        return;
    }

    double FakeTrackSignalGenerator::GetStartPitchMax() const
    {
        return fStartPitchMax;
    }

    void FakeTrackSignalGenerator::SetStartPitchMax( double aPitchMax )
    {
        fStartPitchMax = aPitchMax;
        return;
    }

    double FakeTrackSignalGenerator::GetStartPitchMin() const
    {
        return fStartPitchMin;
    }

    void FakeTrackSignalGenerator::SetStartPitchMin( double aPitchMin )
    {
        fStartPitchMin = aPitchMin;
        return;
    }

    double FakeTrackSignalGenerator::GetPitchMin() const
    {
        return fPitchMin;
    }

    void FakeTrackSignalGenerator::SetPitchMin( double aPitchMin )
    {
        fPitchMin = aPitchMin;
        return;
    }

    double FakeTrackSignalGenerator::GetStartVPhase() const
    {
        return fStartVPhase;
    }

    void FakeTrackSignalGenerator::SetStartVPhase( double aPhase )
    {
        fStartVPhase = aPhase;
        return;
    }

    double FakeTrackSignalGenerator::GetStartTimeMin() const
    {
        return fStartTimeMin;
    }

    void FakeTrackSignalGenerator::SetStartTimeMin( double aTimeMin )
    {
        fStartTimeMin = aTimeMin;
        return;
    }

    double FakeTrackSignalGenerator::GetStartTimeMax() const
    {
        return fStartTimeMax;
    }

    void FakeTrackSignalGenerator::SetStartTimeMax( double aTimeMax )
    {
        fStartTimeMax = aTimeMax;
        return;
    }

    double FakeTrackSignalGenerator::GetFrequency() const
    {
        return fLO_frequency;
    }

    void FakeTrackSignalGenerator::SetFrequency( double aFrequency )
    {
        fLO_frequency = aFrequency;
        return;
    }

    double FakeTrackSignalGenerator::GetNTracksMean() const
    {
        return fNTracksMean;
    }

    void FakeTrackSignalGenerator::SetNTracksMean( double aNTracksMean )
    {
        fNTracksMean = aNTracksMean;
        return;
    }

    double FakeTrackSignalGenerator::GetBField() const
    {
        return fBField;
    }

    void FakeTrackSignalGenerator::SetBField( double aBField )
    {
        fBField = aBField;
        return;
    }

    double FakeTrackSignalGenerator::GetCoilCurrent() const
    {
        return fCoilCurrent;
    }

    void FakeTrackSignalGenerator::SetCoilCurrent( double aCoilCurrent )
    {
        fCoilCurrent = aCoilCurrent;
        return;
    }

    int FakeTrackSignalGenerator::GetRandomSeed() const
    {
        return fRandomSeed;
    }

    void FakeTrackSignalGenerator::SetRandomSeed( int aRandomSeed )
    {
        fRandomSeed = aRandomSeed;
        return;
    }
    
    double FakeTrackSignalGenerator::GetHydrogenFraction() const
    {
        return fHydrogenFraction;
    }

    void FakeTrackSignalGenerator::SetHydrogenFraction( double aHydrogenFraction )
    {
        fHydrogenFraction = aHydrogenFraction;
        return;
    }

    int FakeTrackSignalGenerator::GetNEvents() const
    {
        return fNEvents;
    }

    void FakeTrackSignalGenerator::SetNEvents( int aNEvents )
    {
        fNEvents = aNEvents;
        return;
    }

    bool FakeTrackSignalGenerator::GetAharmonicCorrection() const
    {
        return fAharmonicCorrection;
    }

    void FakeTrackSignalGenerator::SetAharmonicCorrection( bool aAharmonicCorrection )
    {
        fAharmonicCorrection = aAharmonicCorrection;
        return;
    }

    bool FakeTrackSignalGenerator::GetPitchCorrection() const
    {
        return fPitchCorrection;
    }


    void FakeTrackSignalGenerator::SetPitchCorrection( bool aPitchCorrection )
    {
        fPitchCorrection = aPitchCorrection;
        return;
    }

    bool FakeTrackSignalGenerator::GetSlopeCorrection() const
    {
        return fSlopeCorrection;
    }


    void FakeTrackSignalGenerator::SetSlopeCorrection( bool aSlopeCorrection )
    {
        fSlopeCorrection = aSlopeCorrection;
        return;
    }

    double FakeTrackSignalGenerator::GetTrapLength() const
    {
        return fTrapLength;
    }

    void FakeTrackSignalGenerator::SetTrapLength( double aTrapLength )
    {
        fTrapLength = aTrapLength;
        return;
    }


    Signal::State FakeTrackSignalGenerator::GetDomain() const
    {
        return fRequiredSignalState;
    }

    void FakeTrackSignalGenerator::SetDomain( Signal::State aDomain )
    {
        if( aDomain == fRequiredSignalState ) return;
        fRequiredSignalState = aDomain;  // pls changed == to =.
        if( fRequiredSignalState == Signal::kTime )
        {
            fDoGenerateFunc = &FakeTrackSignalGenerator::DoGenerateTime;
        }
        else if( fRequiredSignalState == Signal::kFreq )
        {
            fDoGenerateFunc = &FakeTrackSignalGenerator::DoGenerateFreq;
        }
        else
        {
            LWARN( lmclog, "Unknown domain requested: " << aDomain );
        }
        return;
    }


    bool FakeTrackSignalGenerator::DoGenerate( Signal* aSignal )
    {
        return (this->*fDoGenerateFunc)( aSignal );
    }


    void FakeTrackSignalGenerator::ReadFile(std::string filename, std::vector<std::pair<double,double> > &data)
    {
        std::ifstream input( filename );
        std::vector<std::pair<double, double> > readData; // energies/ oscillator strengths
        double bufferE, bufferOsc;
        for( std::string line; getline( input, line ); )
        {
            if(line.empty() || line[0] == std::string("#")) continue;
            std::stringstream ss(line);
            ss >> bufferE;
            ss >> bufferOsc;
            readData.push_back(std::make_pair(bufferE, fabs(bufferOsc)));
        }

        //sort data by energy
        std::sort(readData.begin(), readData.end(), [](const std::pair<double,double> & a, const std::pair<double, double> & b) -> bool { return a.first < b.first; });

        //remove duplicates
        auto equal_lambda = [](const std::pair<double, double> &a, const std::pair<double, double> & b) { return a.first == b.first;};
        readData.erase( std::unique( readData.begin(), readData.end(), equal_lambda ), readData.end() );
        data = readData;
    }

    double FakeTrackSignalGenerator::EnergyLossSpectrum(double eLoss, double oscillator_strength)
    {
        double T = rel_energy(fLO_frequency + 50e6, fBField);
        return (LMCConst::E_Rydberg() / eLoss) * oscillator_strength * log(4. * T * eLoss / pow(LMCConst::E_Rydberg(), 3.) ); // Produces energy loss spectrum (N. Buzinsky report Eqn XXX) 
        // NOTE: because this formula depends only on log T, I do NOT update with each change in kinetic energy (ie. from radiative losses). Including these changes may be better

    }

    double FakeTrackSignalGenerator::GetScatteredPitchAngle(double thetaScatter, double pitchAngle, double phi)
    {
        LMCThreeVector v(cos(phi) * sin(thetaScatter), sin(phi) * sin(thetaScatter), cos(thetaScatter)); //random unit vector at angle thetaScatter around z axis
        LMCThreeVector xHat(1.,0.,0.);

        //Rodrigues' formula, rotate into pitch angle
        LMCThreeVector v_new = v * cos(pitchAngle)  + xHat.Cross(v) * sin(pitchAngle)  + xHat * xHat.Dot(v) * (1. - cos(pitchAngle));

        return acos(v_new.Z());
    }

    double FakeTrackSignalGenerator::GetBField(double z) //magnetic field profile (harmonic model)
    {
        return fBField*(1. + pow(z / fTrapLength, 2.));
    }

    double FakeTrackSignalGenerator::GetPitchAngleZ(double theta_i, double B_i, double B_f)
    {
        double sinTheta_f = sin(theta_i) * sqrt(B_f / B_i);
        return asin(sinTheta_f);
    }

    void FakeTrackSignalGenerator::SetInterpolator(gsl_spline*& interpolant, std::vector< std::pair<double, double> > data)
    {
        // Reads in oscillator strength data, fills boost interpolator with corresponding inverse CDF for subsequent inversion sampling
        std::vector<double> energies, oscillator_strengths, energy_loss;
        double fOsc;

        for (auto it = std::make_move_iterator(data.begin()), end = std::make_move_iterator(data.end()); it != end; ++it)
        {
            fOsc = abs(std::move(it->second));
            if(fOsc)
            {
                energies.push_back(std::move(it->first));
                oscillator_strengths.push_back(fOsc);
                energy_loss.push_back(EnergyLossSpectrum(energies.back(), oscillator_strengths.back()));
            }
        }

        std::vector<double> cdf(energy_loss.size());

        for(int i=1;i<cdf.size();++i) //manual trapezoidal rule for cdf integral
            cdf[i] = cdf[i-1] + (energy_loss[i-1] + energy_loss[i]) * (energies[i] - energies[i-1]) / 2.;

        double cdf_end = cdf.back();
        std::transform(cdf.begin(), cdf.end(), cdf.begin(), [cdf_end](double& c){return c/cdf_end;});

        interpolant = gsl_spline_alloc(gsl_interp_cspline, cdf.size());
        gsl_spline_init(interpolant, cdf.data(), energies.data(), energies.size());
    }

    //extrapolate oscillator strength to +1000 eV energy loss using Aseev energy loss formula
    void FakeTrackSignalGenerator::ExtrapolateData(std::vector< std::pair<double, double> > &data, std::array<double, 3> fitPars)
    {
        //rename for convenience
        double A = fitPars[0];
        double w = fitPars[1];
        double e = fitPars[2];
        const double eMax = 1000.;
        const double dE = 1.;
        double oscillator_strength;

        double eBack = data.back().first;
        while( eBack < eMax)
        {
            eBack +=dE;
            oscillator_strength = A * pow(w,2.) / (pow(w,2.) + 4. * pow(e - eBack, 2.));
            data.push_back(std::pair<double, double>(eBack, oscillator_strength));
        }

    }

    double FakeTrackSignalGenerator::rel_cyc(double energy, double b_field) const
    {
        double cyc_freq = LMCConst::Q() * b_field / LMCConst::M_el_kg();
        return cyc_freq / ( 1. + (energy/LMCConst::M_el_eV()) ) / (2. * LMCConst::Pi()); // takes energy in eV, magnetic field in T, returns in Hz
    }


    double FakeTrackSignalGenerator::rel_energy(double frequency, double b_field) const
    {
        double gamma = LMCConst::Q() * b_field / (2. * LMCConst::Pi() * frequency * LMCConst::M_el_kg()); //omega = q B / m Gamma
        return (gamma - 1.) *LMCConst::M_el_eV(); // takes frequency in Hz, magnetic field in T, returns in kinetic energy eV
    }

    double FakeTrackSignalGenerator::GetCorrectedFrequency(double frequency, double radius) const
    {
        double correction_factor = 1.;

        if( (fPitch != LMCConst::Pi() / 2.) && fAharmonicCorrection)
	{
		correction_factor = fAharmonicCorrectionFactor;
	}

	else if ( (fPitch !=  LMCConst::Pi() / 2.) && ( fPitchCorrection == 1 ) )
	{
            correction_factor += 1. / (2. * pow(tan(fPitch), 2.));

            correction_factor -= 1. / 2. * pow(radius / fTrapLength, 2.);
	}

        return frequency * correction_factor;

    } // non-90 pitches change average B, changing apparent frequency (A. Astari Esfahani et al. (2019) (Eqn 56))

    double FakeTrackSignalGenerator::WaveguidePowerCoupling(double frequency, double pitchAngle)
    {
        const double tWaveguideRadius = 0.00502920;
        double f_c = 1.8412 * LMCConst::C() / (2 * LMCConst::Pi() * tWaveguideRadius);
        double v_phase = LMCConst::C() / sqrt( 1 - pow(f_c / frequency, 2));
        double k_lambda = 2. * LMCConst::Pi() * frequency / v_phase;
        double zMax = 0;
        if (pitchAngle != LMCConst::Pi() / 2.)
            zMax = fTrapLength / tan(pitchAngle);
        return j0(k_lambda * zMax);
    }

    double FakeTrackSignalGenerator::Z11(double aFrequency, double aRadius)
    {
        const double tWaveguideRadius = 0.00502920;
        double f_c = 1.8412 * LMCConst::C() / (2. * LMCConst::Pi() * tWaveguideRadius);

        double tEnergy = rel_energy(aFrequency,fBField-GetTrapField(0,aRadius));
        double tGamma = 1. + tEnergy / LMCConst::M_el_eV();
        double tBeta = sqrt(1. - 1./pow(tGamma,2.));

        double tZ11 = pow(tBeta,2.) / sqrt(1. - pow(f_c / aFrequency,2.));
        return tZ11;
    }

    double FakeTrackSignalGenerator::RadialPowerCoupling(double radius)
    {
        const double tWaveguideRadius = 0.00502920;
        double k_c = 1.8412 / tWaveguideRadius;
        double x = k_c * radius;

        double corr2 = 1. / 4. * pow(j0(x) - jn(2, x), 2.);
        if(x > 1e-100)
            corr2 += pow(j1(x) / x , 2.);
        else
            corr2 += 1./4.; //avoid nans

        return sqrt(2. * corr2);
    }


    double FakeTrackSignalGenerator::GetEnergyLoss(double u, bool hydrogenScatter)
    {
        //uses interpolated cdf to return random energy
        bool krScatter = !hydrogenScatter;
        int gas_index = int(krScatter);

        return gsl_spline_eval(fInterpolators[gas_index], u, fAccelerators[gas_index]);
    }

    double FakeTrackSignalGenerator::GetAxialFrequency() //angular frequency of axial motion Ali et al. (54) (harmonic trap)
    {
        double gamma = 1. + rel_energy(fStartFrequency, fBField) / LMCConst::M_el_eV();
        double v0 = LMCConst::C() * sqrt( 1. - 1. / pow(gamma, 2.));
        return v0 * sin(fPitch) / fTrapLength;

    }

    double FakeTrackSignalGenerator::GetKa2(double eLoss, double T)  //Generate random momentum transfer (Ka_0)^2. Note rndm generator is encapsulated within function
    {
        double ka2Min = pow(eLoss, 2.) / (4.* LMCConst::E_Rydberg() * T) * (1. + eLoss / (2. *  LMCConst::E_Rydberg()));
        double ka2Max = 4. * T / LMCConst::E_Rydberg()  * (1. - eLoss / (2. * T));
        std::uniform_real_distribution<double> uniform_dist(log(ka2Min), log(ka2Max));
        double lnka2Transfer = uniform_dist(fRandomEngine);
        return exp(lnka2Transfer);

    }

    double FakeTrackSignalGenerator::GetCoilField(double aR0, double aZ0, double aRadius, double aZ) const
    {
        double alpha  = sqrt( pow(aR0,2.) + pow(aRadius,2.) + pow(aZ-aZ0,2.) - 2*aR0*aRadius);
        double beta  = sqrt( pow(aR0,2.) + pow(aRadius,2.) + pow(aZ-aZ0,2.)+ 2*aR0*aRadius);
        double k = sqrt( 1 - pow(alpha/beta, 2.) );

        double tBField = LMCConst::MuNull() * fCoilCurrent / (2. * LMCConst::Pi() * pow(alpha,2.) * beta);
        tBField *= ( ( pow(aR0,2.) - pow(aRadius,2.) - pow(aZ-aZ0,2.) ) * gsl_sf_ellint_Ecomp(k,GSL_PREC_DOUBLE) + pow(alpha,2.) * gsl_sf_ellint_Kcomp(k,GSL_PREC_DOUBLE) );
        
        return tBField;
    }

    void FakeTrackSignalGenerator::AddConst(std::vector<double> &aVec, const double aFactor, const double aConst)
    {
        std::transform(aVec.begin(), aVec.end(), aVec.begin(), [=](double &x) {return aFactor*x+aConst;});
    }

    //Converts quarter-period motion into full period (Aharmonic)
    std::pair<std::vector<double>, std::vector<double> > FakeTrackSignalGenerator::GetFullCycle(std::pair< std::vector<double>, std::vector<double> > aHarmonicSolver)
    {
        std::vector<double> tTimesOrig = aHarmonicSolver.first;
        std::vector<double> tZPositionsOrig = aHarmonicSolver.second;

        std::vector<double> tTimes = tTimesOrig;
        std::vector<double> tZPositions = tZPositionsOrig;
        std::vector<double> tTimesRev = tTimesOrig;
        std::vector<double> tZPositionsRev = tZPositionsOrig;

        tTimesRev.pop_back();
        tZPositionsRev.pop_back();

        std::reverse(tTimesRev.begin(), tTimesRev.end());
        std::reverse(tZPositionsRev.begin(), tZPositionsRev.end());

        AddConst(tTimesRev, -1, 2*tTimes.back());

        //Append electron returning to origin
        tZPositions.insert(tZPositions.end(), tZPositionsRev.begin(), tZPositionsRev.end());
        tTimes.insert(tTimes.end(), tTimesRev.begin(), tTimesRev.end());
        tTimes.pop_back();
        tZPositions.pop_back();

        //Double the whole thing so the e- repeats path over negative z
        std::vector<double> tTimesNeg = tTimes;
        std::vector<double> tZPositionsNeg = tZPositions;

        std::transform(tZPositionsNeg.begin(), tZPositionsNeg.end(), tZPositionsNeg.begin(), [](double &c) {return -c;});
        AddConst(tTimesNeg, 1, tTimesNeg.back());

        tZPositions.insert(tZPositions.end(), tZPositionsNeg.begin(), tZPositionsNeg.end());
        tTimes.insert(tTimes.end(), tTimesNeg.begin(), tTimesNeg.end());

        return std::pair< std::vector<double>, std::vector<double> >(tTimes, tZPositions);
    }

    double FakeTrackSignalGenerator::AharmonicPowerCoupling(double aRadius, double aTheta)
    {
        auto tAharmonicSolver = GetParticleTimes(aRadius, aTheta);
        tAharmonicSolver = GetFullCycle(tAharmonicSolver);
        std::vector<double> tTimes = tAharmonicSolver.first;
        std::vector<double> tZPositions = tAharmonicSolver.second;
        std::vector<double> tBFields;

        for(unsigned i=0;i<tZPositions.size();++i)
            tBFields.push_back(fBField - GetTrapField(tZPositions[i],aRadius));

        double tBAverage = fBField - GetAverageMagneticField(aRadius, aTheta);
        double tEnergy = rel_energy(fStartFrequency,fBField-GetTrapField(0,aRadius));
        double tGamma = 1. + tEnergy / LMCConst::M_el_eV();

        double tRatioBToFreq = LMCConst::Q() / (LMCConst::M_el_kg() * tGamma);
        double tOmega0 = tBAverage * tRatioBToFreq;

        std::vector<double> tOmegas;
        for(unsigned i=0;i<tZPositions.size();++i)
            tOmegas.push_back(tBFields[i] * tRatioBToFreq );

        const double tWaveguideRadius = 0.00502920;
        double f_c = 1.8412 * LMCConst::C() / (2 * LMCConst::Pi() * tWaveguideRadius);
        double v_phase = LMCConst::C() / sqrt( 1 - pow(f_c / fStartFrequency, 2));
        double k_lambda = 2. * LMCConst::Pi() * fStartFrequency / v_phase;

        std::vector<double> tTimeDiffs;
        for(unsigned i=0; i<tTimes.size() - 1;++i)
            tTimeDiffs.push_back(tTimes[i+1] - tTimes[i]);


        std::vector<double> tPhi={0};
        double tPhiSum = 0.;
        for(unsigned i=0;i<tOmegas.size()-1;++i)
        {
            tPhiSum += tOmegas[i] * tTimeDiffs[i];
            tPhi.push_back(tPhiSum);
        }

        double tCosSum = 0.;
        double tSinSum = 0.;

        for(unsigned i=0; i<tTimes.size()-1;++i)
        {
            tCosSum += cos( tPhi[i] + k_lambda * tZPositions[i]  - tOmega0 * tTimes[i] ) * tTimeDiffs[i];
            tSinSum += sin( tPhi[i] + k_lambda * tZPositions[i]  - tOmega0 * tTimes[i] ) * tTimeDiffs[i];

        }

        double tA0 = sqrt(pow(tCosSum,2.) + pow(tSinSum,2.))/ tTimes.back();
        return tA0;
    }

    double FakeTrackSignalGenerator::GetTrapField(double aZ, double aRadius) const
    {
        const double R = 0.01486 / 2.; //Inner Radius of the first round of coils
        const double OD = 0.0004572; //Wires Outer Radius
        const double L = 0.00762; //Length of the Coils
        const double a = 1.006e-2/2; //Waveguide Radius

        std::vector<double> tCoilSeriesRadius = {R + OD / 2., R + OD / 2. * ( 1 + sqrt(3)), R + OD / 2. * ( 1. + 2. * sqrt(3.)), R + OD / 2. * ( 1. + 3.*sqrt(3.))};

        const unsigned nPositions = 16;
        std::vector<double> tCoilSeriesPositions(nPositions);

        for(unsigned i=0; i<nPositions; ++i)
            tCoilSeriesPositions[i] = - L/2. + (2*i+1) * OD/2.;

        double B_value = 0;

        for( unsigned i=0; i<tCoilSeriesRadius.size(); ++i)
        {
            for( unsigned j=0; j<tCoilSeriesPositions.size(); ++j)
            {
                B_value += GetCoilField(tCoilSeriesRadius[i], tCoilSeriesPositions[j], aRadius, aZ);
            }
        }

        return B_value;
    }

    double FakeTrackSignalGenerator::GetZMax(double aTheta, double aRadius) const
    {
        double tBMax =  (fBField - GetTrapField(0,aRadius)) / pow(sin(aTheta), 2.);

        const double tZMaxMax = 0.02; //max possible z max
        const unsigned tNPoints = 10000;
        const double dZ = tZMaxMax / tNPoints;

        std::vector<double> tBDifference;
        for(double z=0; z < tZMaxMax; z+=dZ)
            tBDifference.push_back(fabs( -tBMax + fBField - GetTrapField(z,aRadius)) );

        int tTurningPointIndex = std::min_element(tBDifference.begin(),tBDifference.end()) - tBDifference.begin();

        return tTurningPointIndex * dZ;
    }

    std::pair< std::vector<double>, std::vector<double> > FakeTrackSignalGenerator::GetParticleTimes(double aRadius, double aTheta) const
    {
        double tTime = 0.;
        std::vector<double> tTimes, tZPositions;

        const int nPoints = 2000;
        double tZMax = GetZMax(aTheta, aRadius);
        const double dZ = tZMax / nPoints;
        double xDummy;
        double tDummy = 0;
        double tEnergy = rel_energy(fStartFrequency,fBField-GetTrapField(0,aRadius));
        double tGamma = 1. + tEnergy / LMCConst::M_el_eV();
        double tBeta = sqrt(1. - 1./pow(tGamma,2.));
        double v0 = tBeta * LMCConst::C();

        for(unsigned i=0; i<nPoints; ++i)
        {
            double tZi = i * dZ; 

            xDummy = v0 * sqrt(1. - pow(sin(aTheta),2.) * (fBField - GetTrapField(tZi, aRadius)) / (fBField - GetTrapField(0,aRadius)));
            tDummy += dZ / xDummy;
            if(isnan(tDummy)) break;

            tZPositions.push_back(tZi);
            tTimes.push_back(tDummy);
        }

        return std::pair< std::vector<double>, std::vector<double> >(tTimes, tZPositions);
    }

    double FakeTrackSignalGenerator::GetAverageMagneticField(double aRadius, double aTheta) const
    {

        auto tAharmonicSolver = GetParticleTimes(aRadius, aTheta);
        std::vector<double> tTimes = tAharmonicSolver.first;
        std::vector<double> tZPositions = tAharmonicSolver.second;
        std::vector<double> tBFields;

        for(unsigned i=0;i<tZPositions.size();++i)
            tBFields.push_back(GetTrapField(tZPositions[i],aRadius));

        double tTotalSum = 0.;

        for(unsigned i=0; i<tTimes.size() - 1; ++i)
        {
            tTotalSum +=  tBFields[i] * ( tTimes[i+1] - tTimes[i]);
        }

       return tTotalSum / tTimes.back();
    }

    
    void FakeTrackSignalGenerator::SetTrackProperties(Track &aTrack, int TrackID, double aTimeOffset)
    {
        double current_energy = 0.;
        double energy_loss = 0.;
        double new_energy = 0.;
        double scattering_cdf_val = 0.;
        double zScatter = 0.;
        bool scatter_hydrogen;
        double pitch_angle;
        double theta_scatter;
        const double deg_to_rad = LMCConst::Pi() / 180.;

        std::uniform_real_distribution<double> starttime_distribution(fStartTimeMin,fStartTimeMax);
        std::uniform_real_distribution<double> startpitch_distribution(cos(fStartPitchMin * deg_to_rad),cos(fStartPitchMax * deg_to_rad));
        std::uniform_real_distribution<double> dist(0,1);

        if(TrackID==0)
        {
            if (aTimeOffset==0) // first event
        	{
                fStartTime = starttime_distribution(fRandomEngine);
        	}
            else
            {
                fStartTime = aTimeOffset;
            }

            if(fUseFrequencyDistribution)
            {
                fStartFrequency = fStartFrequencyDistribution->Generate();
            }
            else
            {
                fStartFrequency = rel_cyc(fStartEnergyDistribution->Generate(), fBField);
            }

            do
            {
                fPitch = acos(startpitch_distribution(fRandomEngine));
                double z0 = fz0Distribution->Generate();
                fPitch = GetPitchAngleZ(fPitch, GetBField(z0), fBField);

            } while(fPitch < fPitchMin * deg_to_rad );

            fRadius = fRadiusDistribution->Generate();

            if(fAharmonicCorrection)
            {
                double tBDifference = GetAverageMagneticField(fRadius, fPitch) - GetTrapField(0,0);
                fAharmonicCorrectionFactor = 1. - tBDifference / fBField;
            }
            aTrack.Radius = fRadius;
            aTrack.StartTime = fStartTime;
            aTrack.StartFrequency = fStartFrequency;
        }

       else
       {
            fStartTime = fEndTime + 0.;  // old track endtime + margin=0.
            current_energy = rel_energy(fCurrentFrequency, fBField); // convert current frequency to energy

            scatter_hydrogen = ( dist(fRandomEngine) <= fHydrogenFraction); // whether to scatter of H2 in this case
            scattering_cdf_val = dist(fRandomEngine); // random continous variable for scattering inverse cdf input
            energy_loss = GetEnergyLoss(scattering_cdf_val, scatter_hydrogen); // get a random energy loss using the inverse sampling theorem, scale to eV
            theta_scatter = fScatteringAngleDistribution->Generate(); // get scattering angle (NOT Pitch)

            // Compute new pitch angle, given initial pitch angle, scattering angle. Account for how kinematics change with different axial position of scatter
            if(fPitch != LMCConst::Pi() / 2.)
                zScatter = fTrapLength / tan(fPitch) * sin( GetAxialFrequency() * fStartTime);

            double thetaTop = GetPitchAngleZ(fPitch, fBField, GetBField(zScatter));
            double newThetaTop = GetScatteredPitchAngle( theta_scatter, thetaTop, 2. * LMCConst::Pi() * dist(fRandomEngine) ); //Get pitch angle
            fPitch = GetPitchAngleZ(newThetaTop, GetBField(zScatter), fBField);

            new_energy = current_energy - energy_loss; // new energy after loss, in eV
            fStartFrequency = rel_cyc(new_energy, fBField);
            fCurrentFrequency = fStartFrequency;
            if(fAharmonicCorrection)
            {
                double tBDifference = GetAverageMagneticField(fRadius, fPitch) - GetTrapField(0,0);
                fAharmonicCorrectionFactor = 1. - tBDifference / fBField;
            }

            aTrack.StartTime = fEndTime + 0.; // margin of time is 0.
            aTrack.StartFrequency = fStartFrequency;
        }

        fSlope = fSlopeDistribution->Generate();
        if(fSlopeCorrection) fSlope *= pow(RadialPowerCoupling(fRadius),2.);

        double coupling;
        if(fAharmonicCorrection)
        {
            double tZRatio = Z11(fStartFrequency,fRadius)/Z11(rel_cyc(18600., fBField),0);
            coupling = AharmonicPowerCoupling(fRadius, fPitch) * RadialPowerCoupling(fRadius) * sqrt(tZRatio);
            fAharmonicPowerCoupling = coupling;

        }
        else
        {
            coupling = WaveguidePowerCoupling(fStartFrequency, fPitch) * RadialPowerCoupling(fRadius);
        }

        fTrackLength = fTrackLengthDistribution->Generate();
        fEndTime = fStartTime + fTrackLength;  // reset endtime.
        aTrack.Slope = fSlope;
        aTrack.TrackLength = fTrackLength;
        aTrack.EndTime = aTrack.StartTime + aTrack.TrackLength;
        aTrack.LOFrequency = fLO_frequency;
        aTrack.TrackPower = fSignalPower * pow(coupling,2.);
        aTrack.StartFrequency = GetCorrectedFrequency(aTrack.StartFrequency, aTrack.Radius);
        aTrack.PitchAngle = fPitch * 180. / LMCConst::Pi();
    }

    void FakeTrackSignalGenerator::InitiateEvent(Event* anEvent, int eventID)
    {
        int random_seed_val;
        if ( fRandomSeed != 0 )
        {
            random_seed_val = fRandomSeed;
        }
        else
        {
            std::random_device rd;
            random_seed_val = rd();
        }

        std::geometric_distribution<int> ntracks_distribution(1./fNTracksMean);
        fNTracks = ntracks_distribution(fRandomEngine)+1;

        anEvent->fEventID = eventID;
        anEvent->fLOFrequency = fLO_frequency;
        anEvent->fRandomSeed = random_seed_val;
    }


    bool FakeTrackSignalGenerator::DoGenerateTime( Signal* aSignal )
    {

    	FileWriter* aRootTreeWriter = RootTreeWriter::get_instance();
    	aRootTreeWriter->SetFilename(fRootFilename);
        aRootTreeWriter->OpenFile("RECREATE");

        const unsigned nChannels = fNChannels;
        const double tLocustStep = 1./aSignal->DecimationFactor()/(fAcquisitionRate*1.e6);
        double signalAmplitude;
        double tTimeOffset = 0;

        for(unsigned eventID = 0; eventID < fNEvents; ++eventID)// event loop.
        {
            Event* anEvent = new Event();
            InitiateEvent(anEvent, eventID);
            Track aTrack;
            SetTrackProperties(aTrack, 0, tTimeOffset);
            anEvent->AddTrack(aTrack);

            unsigned tTrackIndex = 0;

            while( tTrackIndex < anEvent->fNTracks) //loop over tracks in event
            {

                for(unsigned ch = 0; ch < nChannels; ++ch) // over all channels
                {
                	double LO_phase = 0.;
                    double voltage_phase = fStartVPhase;
                    unsigned tTrackIndexRange[2] = {static_cast<unsigned>(fStartTime / tLocustStep), static_cast<unsigned>(fEndTime / tLocustStep)};
                    tTrackIndexRange[1] = std::min(tTrackIndexRange[1], aSignal->TimeSize()*aSignal->DecimationFactor());
                    fCurrentFrequency = fStartFrequency;

                    for( unsigned index = tTrackIndexRange[0]; index < tTrackIndexRange[1]; ++index ) // advance sampling time
                    {
                        LO_phase += 2.*LMCConst::Pi()*fLO_frequency * tLocustStep;
                        fCurrentFrequency += fSlope * 1.e6/1.e-3 * tLocustStep;
                        voltage_phase += 2.*LMCConst::Pi()*GetCorrectedFrequency(fCurrentFrequency, aTrack.Radius) * tLocustStep;
                        if(fAharmonicCorrection)
                        {
                            signalAmplitude = sqrt(50.) * sqrt(fSignalPower) * fAharmonicPowerCoupling;
                        }
                        else
                        {
                            signalAmplitude = sqrt(50.) * sqrt(fSignalPower) * WaveguidePowerCoupling(fCurrentFrequency, fPitch) * RadialPowerCoupling(fRadius);
                        }
                        aSignal->LongSignalTimeComplex()[ch*aSignal->TimeSize()*aSignal->DecimationFactor() + index][0] += signalAmplitude * cos(voltage_phase-LO_phase);
                        aSignal->LongSignalTimeComplex()[ch*aSignal->TimeSize()*aSignal->DecimationFactor() + index][1] += signalAmplitude * cos(-LMCConst::Pi()/2. + voltage_phase-LO_phase);
                    }

                } //channel loop

                ++tTrackIndex;
                tTimeOffset = fEndTime;
                SetTrackProperties(aTrack, tTrackIndex, tTimeOffset);
                double tPitchMinRad = fPitchMin * LMCConst::Pi() / 180.;

                if( (!fNTracksMean && (fPitch < tPitchMinRad )) || (fNTracksMean && (tTrackIndex == fNTracks)))
                {
                    break;
                }
                else
                {
                    anEvent->AddTrack(aTrack);
                }

            } //track loop

            aRootTreeWriter->WriteEvent(anEvent);
            delete anEvent;
        } //event loop

        aRootTreeWriter->CloseFile();
        return true;
    }

    bool FakeTrackSignalGenerator::DoGenerateFreq( Signal* aSignal )
    {
        return true;
    }

} /* namespace locust */
