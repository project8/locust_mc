/*
 * LMCTestFIRFilter.cc
 *
 *  Created on: May 3 2019
 *      Author: plslocum
 */

#include <cmath>
#include "LMCTestFIRFilterGenerator.hh"
#include "LMCDigitizer.hh"
#include "logger.hh"
#include "LMCConst.hh"


using std::string;

namespace locust
{
    LOGGER( lmclog, "TestFIRFilterGenerator" );

    MT_REGISTER_GENERATOR(TestFIRFilterGenerator, "test-firfilter");

    TestFIRFilterGenerator::TestFIRFilterGenerator( const std::string& aName ) :
        Generator( aName ),
        fDoGenerateFunc( &TestFIRFilterGenerator::DoGenerateTime ),
        gfilter_filename("blank.txt"),
        fFilter_resolution( 0. ),
        fLO_frequency( 20.05e9 ),
        fRF_frequency( 20.1e9 ),
        fAmplitude( 5.e-8 )
    {
        fRequiredSignalState = Signal::kTime;
    }

    TestFIRFilterGenerator::~TestFIRFilterGenerator()
    {
    }


    bool TestFIRFilterGenerator::Configure( const scarab::param_node* aParam )
    {
        if( aParam == NULL) return true;

        if( aParam->has( "filter-filename" ) )
        {
            gfilter_filename = aParam->get_value< std::string >( "filter-filename" );
        }

        if( aParam->has( "filter-resolution" ) )
        {
            fFilter_resolution = aParam->get_value< double >( "filter-resolution" );
        }


        if( aParam->has( "rf-frequency" ) )
        {
        SetRFFrequency( aParam->get_value< double >( "rf-frequency", fRF_frequency ) );
        }

        if( aParam->has( "lo-frequency" ) )
        {
        SetLOFrequency( aParam->get_value< double >( "lo-frequency", fLO_frequency ) );
        }


        if( aParam->has( "amplitude" ) )
        {
        SetAmplitude( aParam->get_value< double >( "amplitude", fAmplitude ) );
        }


        if( aParam->has( "domain" ) )
        {
            string domain = aParam->get_value( "domain" );
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


        return true;



    }


    void TestFIRFilterGenerator::Accept( GeneratorVisitor* aVisitor ) const
    {
        aVisitor->Visit( this );
        return;
    }

    double TestFIRFilterGenerator::GetRFFrequency() const
    {
        return fRF_frequency;
    }

    void TestFIRFilterGenerator::SetRFFrequency( double aFrequency )
    {
        fRF_frequency = aFrequency;
        return;
    }

    double TestFIRFilterGenerator::GetLOFrequency() const
    {
        return fLO_frequency;
    }

    void TestFIRFilterGenerator::SetLOFrequency( double aFrequency )
    {
        fLO_frequency = aFrequency;
        return;
    }

    double TestFIRFilterGenerator::GetAmplitude() const
    {
        return fAmplitude;
    }

    void TestFIRFilterGenerator::SetAmplitude( double aAmplitude )
    {
        fAmplitude = aAmplitude;
        return;
    }



    Signal::State TestFIRFilterGenerator::GetDomain() const
    {
        return fRequiredSignalState;
    }

    void TestFIRFilterGenerator::SetDomain( Signal::State aDomain )
    {
        if( aDomain == fRequiredSignalState ) return;
        fRequiredSignalState = aDomain;  // pls changed == to =.
        if( fRequiredSignalState == Signal::kTime )
        {
            fDoGenerateFunc = &TestFIRFilterGenerator::DoGenerateTime;
        }
        else if( fRequiredSignalState == Signal::kFreq )
        {
            fDoGenerateFunc = &TestFIRFilterGenerator::DoGenerateFreq;
        }
        else
        {
            LWARN( lmclog, "Unknown domain requested: " << aDomain );
        }
        return;
    }


    double* TestFIRFilterGenerator::GetFIRFilter(int nskips)
    {

    FILE *fp;
    double *filterarray = new double[1000];
    double filter;
    double index;
    fp = fopen(gfilter_filename.c_str(),"r");
    int count = 0;


    for (int i=0; i<1000; i++)
      filterarray[i] = -99.;



      while (!feof(fp))
        {
        fscanf(fp, "%lf %lf\n", &index, &filter);
        if (count%nskips==0) filterarray[count/nskips] = filter;
    //    printf("filter %d is %g\n", count, filterarray[count]);
        count += 1;
        }

    fclose(fp);
    return filterarray;

    }

    int TestFIRFilterGenerator::GetNFilterBins(double* filterarray)
    {
    int nbins = 0;
    for (int i=0; i<1000; i++)
      {
      if (filterarray[i]>0.) nbins += 1;
      }    
    return nbins;
    }

    double TestFIRFilterGenerator::GetFIRSample(double* filterarray, int nfilterbins, double dtfilter, double fieldamplitude, double fieldphase, double fieldfrequency)
    {
    	std::deque<double> fieldarrived = {0.};

    	for (int i=0; i<nfilterbins; i++)  // initialize field and coefficients.
    	  {
    	  fieldarrived.emplace(fieldarrived.begin()+i, 0.);  // size deque array
    	  }

    	for (int i=0; i<nfilterbins; i++)  // populate filter.
    	  {
    	  fieldphase += 2.*3.1415926*fieldfrequency*dtfilter;
    	  fieldarrived.push_back(fieldamplitude*cos(fieldphase));
    	  fieldarrived.pop_front();
    	  }

    	double total = 0.;
    	for (int j=0; j<nfilterbins; j++)  // sum products in filter.
    	  {
    	  total += fieldarrived.at(j)*filterarray[j];
    	  }

//printf("total is %g\n", total); getchar();
    return total;
    }



    bool TestFIRFilterGenerator::DoGenerate( Signal* aSignal )
    {
        return (this->*fDoGenerateFunc)( aSignal );
    }

    bool TestFIRFilterGenerator::DoGenerateTime( Signal* aSignal )
    {

        const unsigned nchannels = fNChannels;

        double LO_phase = 0.;
        double field_phase = 0.;
        double VoltageSample = 0.;
        double* filterarray = GetFIRFilter(1);
        int nfilterbins = GetNFilterBins(filterarray);
        double dtfilter = fFilter_resolution;

        for (unsigned ch = 0; ch < nchannels; ++ch)
        {
            for( unsigned index = 0; index < aSignal->TimeSize()*aSignal->DecimationFactor(); ++index )
            {

                LO_phase = 2.*LMCConst::Pi()*fLO_frequency*(double)index/aSignal->DecimationFactor()/(fAcquisitionRate*1.e6);
                field_phase = 2.*LMCConst::Pi()*fRF_frequency*(double)index/aSignal->DecimationFactor()/(fAcquisitionRate*1.e6);
                VoltageSample = GetFIRSample(filterarray, nfilterbins, dtfilter, fAmplitude, field_phase, fRF_frequency);

// factor of 2 is needed for cosA*cosB = 1/2*(cos(A+B)+cos(A-B)); usually we leave out the 1/2 for e.g. sinusoidal RF.
// This allows for correct gain in Locust-Katydid analysis chain.
                aSignal->LongSignalTimeComplex()[ch*aSignal->TimeSize()*aSignal->DecimationFactor() + index][0] += 2.*VoltageSample*cos(LO_phase);  
                aSignal->LongSignalTimeComplex()[ch*aSignal->TimeSize()*aSignal->DecimationFactor() + index][1] += 2.*VoltageSample*cos(LMCConst::Pi()/2. + LO_phase);

//printf("signal %d is with acqrate %g, lo %g and rf %g is %g\n", index, fAcquisitionRate, fLO_frequency, fRF_frequency, aSignal->LongSignalTimeComplex()[ch*aSignal->TimeSize()*aSignal->DecimationFactor() + index][0]); getchar();


            }
        }
        return true;
    }

    bool TestFIRFilterGenerator::DoGenerateFreq( Signal* aSignal )
    {
        return true;
    }

} /* namespace locust */
