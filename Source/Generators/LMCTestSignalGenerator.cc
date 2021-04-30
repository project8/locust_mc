/*
 * LMCTestSignal.cc
 *
 *  Created on: Jan 14 2015
 *      Author: plslocum after nsoblath
 */

#include <cmath>
#include "LMCTestSignalGenerator.hh"
#include "LMCDigitizer.hh"
#include "logger.hh"
#include "LMCConst.hh"


using std::string;

namespace locust
{
    LOGGER( lmclog, "TestSignalGenerator" );

    MT_REGISTER_GENERATOR(TestSignalGenerator, "test-signal");

    TestSignalGenerator::TestSignalGenerator( const std::string& aName ) :
        Generator( aName ),
        fDoGenerateFunc( &TestSignalGenerator::DoGenerateTime ),
        fLO_frequency( 20.05e9 ),
        fRF_frequency( 20.1e9 ),
		fAMModfrequency( 0.0e6 ),
		fAMModSelect( 0 ),
        fAmplitude( 5.e-8 ),
		fMixingProduct( false ),
    	fWriteRootHisto( false ),
		fWriteRootGraph( false ),
		fRootFilename( "LocustTestSignal.root")
    {
        fRequiredSignalState = Signal::kTime;
    }

    TestSignalGenerator::~TestSignalGenerator()
    {
    }


    bool TestSignalGenerator::Configure( const scarab::param_node& aParam )
    {
        if( aParam.has( "rf-frequency" ) )
        {
            SetRFFrequency( aParam.get_value< double >( "rf-frequency", fRF_frequency ) );
        }

        if( aParam.has( "lo-frequency" ) )
        {
            SetLOFrequency( aParam.get_value< double >( "lo-frequency", fLO_frequency ) );
        }
        if( aParam.has( "am-frequency" ) )
        {
            SetAMModFrequency( aParam.get_value< double >( "am-frequency", fAMModfrequency ) );
        }
        if( aParam.has( "am-mod-select" ) )
        {
            SetAMModSelect( aParam.get_value< int >( "am-mod-select", fAMModSelect ) );
        }
        if( aParam.has( "amplitude" ) )
        {
            SetAmplitude( aParam.get_value< double >( "amplitude", fAmplitude ) );
        }

        if( aParam.has( "mixing-product" ) )
        {
        	SetMixingProduct( aParam.get_value< bool >( "mixing-product", fMixingProduct ));
        }

        if( aParam.has( "write-root-histo" ) )
        {
        	fWriteRootHisto = aParam.get_value< bool >( "write-root-histo", fWriteRootHisto );
        }

        if( aParam.has( "write-root-graph" ) )
        {
        	fWriteRootGraph = aParam.get_value< bool >( "write-root-graph", fWriteRootGraph );
        }

    	if( aParam.has( "root-filename" ) )
        {
            fRootFilename = aParam["root-filename"]().as_string();
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


        return true;



    }


    void TestSignalGenerator::Accept( GeneratorVisitor* aVisitor ) const
    {
        aVisitor->Visit( this );
        return;
    }

    double TestSignalGenerator::GetRFFrequency() const
    {
        return fRF_frequency;
    }

    void TestSignalGenerator::SetRFFrequency( double aFrequency )
    {
        fRF_frequency = aFrequency;
        return;
    }

    double TestSignalGenerator::GetLOFrequency() const
    {
        return fLO_frequency;
    }

    void TestSignalGenerator::SetLOFrequency( double aFrequency )
    {
        fLO_frequency = aFrequency;
        return;
    }

    double TestSignalGenerator::GetAMModFrequency() const
    {
        return fAMModfrequency;
    }

    void TestSignalGenerator::SetAMModFrequency( double aFrequency )
    {
        fAMModfrequency = aFrequency;
        return;
    }

    int TestSignalGenerator::GetAMModSelect() const
    {
        return fAMModSelect;
    }

    void TestSignalGenerator::SetAMModSelect( int aModSelection )
    {
        fAMModSelect = aModSelection;
        return;
    }

    double TestSignalGenerator::GetAmplitude() const
    {
        return fAmplitude;
    }

    void TestSignalGenerator::SetAmplitude( double aAmplitude )
    {
        fAmplitude = aAmplitude;
        return;
    }

    bool TestSignalGenerator::GetMixingProduct() const
    {
        return fMixingProduct;
    }

    void TestSignalGenerator::SetMixingProduct( bool aMixingProduct )
    {
        fMixingProduct = aMixingProduct;
        return;
    }


    Signal::State TestSignalGenerator::GetDomain() const
    {
        return fRequiredSignalState;
    }

    void TestSignalGenerator::SetDomain( Signal::State aDomain )
    {
        if( aDomain == fRequiredSignalState ) return;
        fRequiredSignalState = aDomain;
        if( fRequiredSignalState == Signal::kTime )
        {
            fDoGenerateFunc = &TestSignalGenerator::DoGenerateTime;
        }
        else if( fRequiredSignalState == Signal::kFreq )
        {
            fDoGenerateFunc = &TestSignalGenerator::DoGenerateFreq;
        }
        else
        {
            LWARN( lmclog, "Unknown domain requested: " << aDomain );
        }
        return;
    }

    bool TestSignalGenerator::WriteRootHisto()
    {
		#ifdef ROOT_FOUND
    	FileWriter* aRootHistoWriter = RootHistoWriter::get_instance();
    	aRootHistoWriter->SetFilename(fRootFilename);
    	aRootHistoWriter->OpenFile("UPDATE");
    	TH1D* aHisto = new TH1D("testhisto", "histotitle", 200, 0., 400.);
    	for (unsigned i=0; i<200; i++)
    	{
    		aHisto->SetBinContent(i+1, (double)i);
    	}
    	aRootHistoWriter->Write1DHisto(aHisto);
        aRootHistoWriter->CloseFile();
		#endif
        return true;
    }

    bool TestSignalGenerator::WriteRootGraph()
    {
		#ifdef ROOT_FOUND
    	FileWriter* aRootGraphWriter = RootGraphWriter::get_instance();
    	aRootGraphWriter->SetFilename(fRootFilename);
    	aRootGraphWriter->OpenFile("UPDATE");
    	std::vector<double> xVector{1,2,3,4,5};
    	std::vector<double> yVector{1,2,3,4,5};
    	aRootGraphWriter->WriteVectorGraph(xVector, yVector);
        aRootGraphWriter->CloseFile();
		#endif
        return true;
    }




    bool TestSignalGenerator::DoGenerate( Signal* aSignal )
    {
        return (this->*fDoGenerateFunc)( aSignal );
    }

    bool TestSignalGenerator::DoGenerateTime( Signal* aSignal )
    {

    	if (fWriteRootHisto) WriteRootHisto();
    	if (fWriteRootGraph) WriteRootGraph();

        const unsigned nchannels = fNChannels;

        double LO_phase = 0.;
        double voltage_phase = 0.;
        double AMMod_phase = 0.;

        for (unsigned ch = 0; ch < nchannels; ++ch)
        {
            for( unsigned index = 0; index < aSignal->TimeSize()*aSignal->DecimationFactor(); ++index )
            {

                LO_phase = 2.*LMCConst::Pi()*fLO_frequency*(double)index/aSignal->DecimationFactor()/(fAcquisitionRate*1.e6);
                voltage_phase = 2.*LMCConst::Pi()*fRF_frequency*(double)index/aSignal->DecimationFactor()/(fAcquisitionRate*1.e6);
                AMMod_phase = 2.*LMCConst::Pi()*fAMModfrequency*(double)index/aSignal->DecimationFactor()/(fAcquisitionRate*1.e6);

                if (fMixingProduct)  // keep upper and lower sidebands RF+LO and RF-LO
                {
        		    aSignal->LongSignalTimeComplex()[ch*aSignal->TimeSize()*aSignal->DecimationFactor() + index][0] +=  2. * sqrt(50.)*fAmplitude*pow(cos(AMMod_phase),fAMModSelect)*cos(voltage_phase) * sin(LO_phase);
        		    aSignal->LongSignalTimeComplex()[ch*aSignal->TimeSize()*aSignal->DecimationFactor() + index][1] +=  2. * sqrt(50.)*fAmplitude*pow(cos(AMMod_phase),fAMModSelect)*cos(voltage_phase) * cos(LO_phase);
                }
                else // keep only lower sideband RF-LO
                {
                	aSignal->LongSignalTimeComplex()[ch*aSignal->TimeSize()*aSignal->DecimationFactor() + index][0] += sqrt(50.)*fAmplitude*pow(cos(AMMod_phase),fAMModSelect)*cos(voltage_phase-LO_phase);
                	aSignal->LongSignalTimeComplex()[ch*aSignal->TimeSize()*aSignal->DecimationFactor() + index][1] += sqrt(50.)*fAmplitude*pow(cos(AMMod_phase),fAMModSelect)*cos(-LMCConst::Pi()/2. + voltage_phase-LO_phase);
                }

/*                if (fMixingProduct)  // keep upper and lower sidebands RF+LO and RF-LO
                {
        		    aSignal->LongSignalTimeComplex()[ch*aSignal->TimeSize()*aSignal->DecimationFactor() + index][0] +=  2. * sqrt(50.)*fAmplitude*pow(cos(carrier_phase),2.) * cos(voltage_phase);
        		    aSignal->LongSignalTimeComplex()[ch*aSignal->TimeSize()*aSignal->DecimationFactor() + index][1] +=  2. * sqrt(50.)*fAmplitude*pow(cos(carrier_phase),2.) * sin(voltage_phase);
                }
                else // keep only lower sideband RF-LO
                {
//                	aSignal->LongSignalTimeComplex()[ch*aSignal->TimeSize()*aSignal->DecimationFactor() + index][0] += sqrt(50.)*fAmplitude*pow(cos(carrier_phase),1.)*cos(voltage_phase-LO_phase);
//                	aSignal->LongSignalTimeComplex()[ch*aSignal->TimeSize()*aSignal->DecimationFactor() + index][1] += sqrt(50.)*fAmplitude*pow(cos(carrier_phase),1.)*cos(-LMCConst::Pi()/2. + voltage_phase-LO_phase);
                }
*/
            }
        }
        return true;
    }

    bool TestSignalGenerator::DoGenerateFreq( Signal* aSignal )
    {
        return true;
    }

} /* namespace locust */
