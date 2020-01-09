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
#include <random>
#include <math.h>
#include <sstream>
#include <string>

using std::string;

namespace locust
{
    LOGGER( lmclog, "FakeTrackSignalGenerator" );

    MT_REGISTER_GENERATOR(FakeTrackSignalGenerator, "fake-track");

    FakeTrackSignalGenerator::FakeTrackSignalGenerator( const std::string& aName ) :
        Generator( aName ),
        fDoGenerateFunc( &FakeTrackSignalGenerator::DoGenerateTime ),
        fSignalPower( 0. ),
        fStartFrequencyMax( 0. ),
        fStartFrequencyMin( 0. ),
        fStartVPhase( 0. ),
        fSlopeMean( 0. ),
        fSlopeStd( 0. ),
        fStartTimeMin( 0. ),
        fStartTimeMax( 0. ),
        fStartPitchMin( 89.9 ),
        fStartPitchMax( 90. ),
        fPitchMin( 0. ),
        fLO_frequency( 0. ),
        fTrackLengthMean( 0. ),
        fNTracksMean( 0. ),
        fBField(1.0),
        fRandomSeed(0),
        fNEvents(1),
        fPitchCorrection( true ),
        fRandomEngine(0),
        fHydrogenFraction(1),
        fTrapLength(0.1784),  //Phase II harmonic trap L0 (A. Ashtari Esfahani et al.- Phys. Rev. C 99, 055501 )
        fH2Interpolant(std::vector<double>(1).data(),std::vector<double>(1).data(),1,0),
        fKrInterpolant(std::vector<double>(1).data(),std::vector<double>(1).data(),1,0),
        fRoot_filename("LocustEvent.root"),
        fSlope( 0. ),
        fPitch( 0. ),
        fTrackLength( 0. ),
        fStartTime( 0. ),
        fEndTime( 0. ),
        fStartFrequency( 0. ),
        fCurrentFrequency( 0. ),
        fNTracks(0)
    {
        fRequiredSignalState = Signal::kTime;
    }

    FakeTrackSignalGenerator::~FakeTrackSignalGenerator()
    {
    }


    bool FakeTrackSignalGenerator::Configure( const scarab::param_node& aParam )
    {
        if( aParam.has( "signal-power" ) )
            SetSignalPower( aParam.get_value< double >( "signal-power", fSignalPower ) );

        if( aParam.has( "start-frequency-max" ) )
            SetStartFrequencyMax( aParam.get_value< double >( "start-frequency-max", fStartFrequencyMax ) );

        if( aParam.has( "start-frequency-min" ) )
            SetStartFrequencyMin( aParam.get_value< double >( "start-frequency-min", fStartFrequencyMin ) );

        if( aParam.has( "start-vphase" ) )
            SetStartVPhase( aParam.get_value< double >( "start-vphase", fStartVPhase ) );

        if( aParam.has( "start-pitch-max" ) )
            SetStartPitchMax( aParam.get_value< double >( "start-pitch-max", fStartPitchMax ) );

        if( aParam.has( "start-pitch-min" ) )
            SetStartPitchMin( aParam.get_value< double >( "start-pitch-min", fStartPitchMin ) );

        if( aParam.has( "min-pitch" ) )
            SetPitchMin( aParam.get_value< double >( "min-pitch", fPitchMin ) );

        if( aParam.has( "slope-mean" ) )
            SetSlopeMean( aParam.get_value< double >( "slope-mean", fSlopeMean ) );

        if( aParam.has( "slope-std" ) )
            SetSlopeStd( aParam.get_value< double >( "slope-std", fSlopeStd ) );

        if( aParam.has( "start-time-max" ) )
            SetStartTimeMax( aParam.get_value< double >( "start-time-max", fStartTimeMax ) );

        if( aParam.has( "start-time-min" ) )
            SetStartTimeMin( aParam.get_value< double >( "start-time-min", fStartTimeMin ) );

        if( aParam.has( "lo-frequency" ) )
            SetFrequency( aParam.get_value< double >( "lo-frequency", fLO_frequency ) );

        if( aParam.has( "track-length-mean" ) )
            SetTrackLengthMean( aParam.get_value< double >( "track-length-mean", fTrackLengthMean ) );

        if (aParam.has( "ntracks-mean") )
            SetNTracksMean( aParam.get_value< double >( "ntracks-mean",fNTracksMean) );

        if (aParam.has("magnetic-field") )
            SetBField(  aParam.get_value< double >("magnetic-field", fBField) );

        if (aParam.has( "random-seed") )
            SetRandomSeed(  aParam.get_value< int >( "random-seed",fRandomSeed) );

        if (aParam.has( "n-events") )
            SetNEvents(  aParam.get_value< int >( "n-events",fNEvents) );

        if (aParam.has( "pitch-correction") )
            SetPitchCorrection(  aParam.get_value< bool >( "pitch-correction", fPitchCorrection) );

        if (aParam.has( "hydrogen-fraction") )
        {
            SetHydrogenFraction(  aParam.get_value< int >( "hydrogen-fraction",fHydrogenFraction) );

            if( fHydrogenFraction > 1 ||  fHydrogenFraction < 0)
                LERROR( lmclog, "hydrogen-fraction must be between 0 and 1!");
        }

        if( aParam.has( "root-filename" ) )
        {
            fRoot_filename = aParam["root-filename"]().as_string();
        }

        if(fRandomSeed)
            fRandomEngine.seed(fRandomSeed);
        else
        {
            std::random_device rd;
            fRandomEngine.seed(rd());
        }

        if(!fNTracksMean && !fPitchMin)
            LERROR( lmclog, "No condition set for NTracks per event! Set one of pitch-min or ntracks-mean");


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
        SetInterpolator(fH2Interpolant,h2Data);
        SetInterpolator(fKrInterpolant,krData);

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

    double FakeTrackSignalGenerator::GetStartFrequencyMax() const
    {
        return fStartFrequencyMax;
    }

    void FakeTrackSignalGenerator::SetStartFrequencyMax( double aFrequencyMax )
    {
        fStartFrequencyMax = aFrequencyMax;
        return;
    }

    double FakeTrackSignalGenerator::GetStartFrequencyMin() const
    {
        return fStartFrequencyMin;
    }

    void FakeTrackSignalGenerator::SetStartFrequencyMin( double aFrequencyMin )
    {
        fStartFrequencyMin = aFrequencyMin;
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

    double FakeTrackSignalGenerator::GetSlopeMean() const
    {
        return fSlopeMean;
    }

    void FakeTrackSignalGenerator::SetSlopeMean( double aSlopeMean )
    {
        fSlopeMean = aSlopeMean;
        return;
    }

    double FakeTrackSignalGenerator::GetSlopeStd() const
    {
        return fSlopeStd;
    }

    void FakeTrackSignalGenerator::SetSlopeStd( double aSlopeStd )
    {
        fSlopeStd = aSlopeStd;
        return;
    }

    double FakeTrackSignalGenerator::GetTrackLengthMean() const
    {
        return fTrackLengthMean;
    }

    void FakeTrackSignalGenerator::SetTrackLengthMean( double aTrackLengthMean )
    {
        fTrackLengthMean = aTrackLengthMean;
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


    bool FakeTrackSignalGenerator::GetPitchCorrection() const
    {
        return fPitchCorrection;
    }


    void FakeTrackSignalGenerator::SetPitchCorrection( bool aPitchCorrection )
    {
        fPitchCorrection = aPitchCorrection;
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
            readData.push_back(std::make_pair(bufferE, bufferOsc));
        }

        //sort data by energy
        std::sort(readData.begin(), readData.end(), [](const std::pair<double,double> & a, const std::pair<double, double> & b) -> bool { return a.first < b.first; });
        data = readData;
    }

    double FakeTrackSignalGenerator::EnergyLossSpectrum(double eLoss, double oscillator_strength)
    {
        double T = rel_energy(fStartFrequencyMax, fBField);
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

    void FakeTrackSignalGenerator::SetInterpolator(boost::math::barycentric_rational<double> &interpolant, std::vector< std::pair<double, double> > data)
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

        interpolant = boost::math::barycentric_rational<double>(cdf.data(), energies.data(), energies.size(), 0); //linear interpolation: don't change, doesnt converge
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

    double FakeTrackSignalGenerator::GetPitchCorrectedFrequency(double frequency) const
    {
        if ( (fPitch !=  LMCConst::Pi() / 2.) && ( fPitchCorrection == 1 ) )
            return frequency * ( 1. + 1. / (2. * pow(tan(fPitch), 2.)));
        else
            return frequency;
    } // non-90 pitches change average B, changing apparent frequency (A. Astari Esfahani et al. (2019) (Eqn 56))

    double FakeTrackSignalGenerator::WaveguidePowerCoupling(double frequency, double pitchAngle)
    {
        double k_lambda = 2. * LMCConst::Pi() * frequency / LMCConst::C();
        double zMax = 0;
        if (pitchAngle != LMCConst::Pi() / 2.)
            zMax = fTrapLength / tan(pitchAngle);
        return j0(k_lambda * zMax);
    }

    double FakeTrackSignalGenerator::GetEnergyLoss(double u, bool hydrogenScatter)
    {
        //uses interpolated cdf to return random energy
        boost::math::barycentric_rational<double> *interpolant;
        if(hydrogenScatter)
        {
            interpolant = &fH2Interpolant;
        }
        else
        {
            interpolant = &fKrInterpolant;
        }
        return (*interpolant)(u);
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
    
    double FakeTrackSignalGenerator::GetThetaScatter(double eLoss, double T)
    {
        double ka2Transfer = GetKa2(eLoss, T);
        double eRatio = eLoss / T;
        double cosTheta = 1. - eRatio /2.  - ka2Transfer  / (2.* (T / LMCConst::E_Rydberg()));
        cosTheta /= sqrt(1. - eRatio);
        return acos(cosTheta);

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

        std::normal_distribution<double> slope_distribution(fSlopeMean,fSlopeStd);
        std::uniform_real_distribution<double> startfreq_distribution(fStartFrequencyMin,fStartFrequencyMax);
        std::exponential_distribution<double> tracklength_distribution(1./fTrackLengthMean);
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
            fStartFrequency = startfreq_distribution(fRandomEngine);
            fPitch = acos(startpitch_distribution(fRandomEngine));
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
            theta_scatter = GetThetaScatter(energy_loss, current_energy); // get scattering angle (NOT Pitch)

            // Compute new pitch angle, given initial pitch angle, scattering angle. Account for how kinematics change with different axial position of scatter
            if(fPitch != LMCConst::Pi() / 2.)
                zScatter = fTrapLength / tan(fPitch) * sin( GetAxialFrequency() * fStartTime);

            double thetaTop = GetPitchAngleZ(fPitch, fBField, GetBField(zScatter));
            double newThetaTop = GetScatteredPitchAngle( theta_scatter, thetaTop, 2. * LMCConst::Pi() * dist(fRandomEngine) ); //Get pitch angle
            fPitch = GetPitchAngleZ(newThetaTop, GetBField(zScatter), fBField);

            new_energy = current_energy - energy_loss; // new energy after loss, in eV
            fStartFrequency = rel_cyc(new_energy, fBField);
            fCurrentFrequency = fStartFrequency;
            aTrack.StartTime = fEndTime + 0.; // margin of time is 0.
            aTrack.StartFrequency = fStartFrequency;
        }

        fSlope = slope_distribution(fRandomEngine);
        fTrackLength = tracklength_distribution(fRandomEngine);
        fEndTime = fStartTime + fTrackLength;  // reset endtime.
        aTrack.Slope = fSlope;
        aTrack.TrackLength = fTrackLength;
        aTrack.EndTime = aTrack.StartTime + aTrack.TrackLength;
        aTrack.LOFrequency = fLO_frequency;
        aTrack.TrackPower = fSignalPower * pow(WaveguidePowerCoupling(fStartFrequency, fPitch),2.);
        aTrack.StartFrequency = GetPitchCorrectedFrequency(aTrack.StartFrequency);
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


    void WriteRootFile(Event* anEvent, TFile* hfile)
    {
    	char buffer[100];
        int n=sprintf(buffer, "Event_%d", anEvent->fEventID);
    	char* treename = buffer;

        TTree *aTree = new TTree(treename,"Locust Tree");
        aTree->Branch("EventID", &anEvent->fEventID, "EventID/I");
        aTree->Branch("ntracks", &anEvent->fNTracks, "ntracks/I");
        aTree->Branch("StartFrequencies", "std::vector<double>", &anEvent->fStartFrequencies);
        aTree->Branch("StartTimes", "std::vector<double>", &anEvent->fStartTimes);
        aTree->Branch("EndTimes", "std::vector<double>", &anEvent->fEndTimes);
        aTree->Branch("TrackLengths", "std::vector<double>", &anEvent->fTrackLengths);
        aTree->Branch("Slopes", "std::vector<double>", &anEvent->fSlopes);
        aTree->Branch("LOFrequency", &anEvent->fLOFrequency, "LOFrequency/D");
        aTree->Branch("RandomSeed", &anEvent->fRandomSeed, "RandomSeed/I");
        aTree->Branch("TrackPower", "std::vector<double>", &anEvent->fTrackPowers);
        aTree->Branch("PitchAngles", "std::vector<double>", &anEvent->fPitchAngles);
        aTree->Fill();
        aTree->Write();
        delete aTree;
    }

    bool FakeTrackSignalGenerator::DoGenerateTime( Signal* aSignal )
    {
        TFile* hfile = new TFile(fRoot_filename.c_str(),"RECREATE");

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
                        voltage_phase += 2.*LMCConst::Pi()*GetPitchCorrectedFrequency(fCurrentFrequency) * tLocustStep;
                        signalAmplitude = sqrt(50.) * sqrt(fSignalPower) * WaveguidePowerCoupling(fCurrentFrequency, fPitch);
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
            WriteRootFile(anEvent, hfile);
            delete anEvent;
        } //event loop

        hfile->Close();
        return true;
    }

    bool FakeTrackSignalGenerator::DoGenerateFreq( Signal* aSignal )
    {
        return true;
    }

} /* namespace locust */
