/*
 * LMCField.cc
 *
 *  Created on: Jun 4, 2021
 *      Author: pslocum
 */

#include "LMCField.hh"


namespace locust
{
    LOGGER( lmclog, "Field" );
    Field::Field():
        fnPixels( 100 ),
        fCentralFrequency(1.63e11),
		fAvgDotProductFactor( 0. ),
        fNModes( 2 ),
        fNChannels( 1 ),
        fbMultiMode( false ),
        fTM111( false ),
        fTE012( false ),
        fTE013( false ),
        fR( 0.18 ),
        fL( 3.0 ),
        fX( 0.010668 ),
        fY( 0.004318 ),
        fCENTER_TO_SHORT( 0.05 ),
        fCENTER_TO_ANTENNA( 0.05 ),
        fPlotModeMaps( false ),
        fOutputPath( TOSTRING(PB_OUTPUT_DIR) )
    {}
    Field::~Field() {}

    FieldCore::FieldCore()
    {}
    FieldCore::~FieldCore() {}



    bool Field::Configure( const scarab::param_node& aParam )
    {

        if( aParam.has( "n-modes" ) )
        {
        	fNModes = aParam["n-modes"]().as_int();
        }

        if( aParam.has( "multi-mode" ) )
        {
    		LPROG(lmclog,"Running in multimode configuration.");
        	fbMultiMode = aParam["multi-mode"]().as_bool();
        }

        if( aParam.has( "tm111-mode" ) )
        {
    		LPROG(lmclog,"Running with TM111 only.");
        	fTM111 = aParam["tm111-mode"]().as_bool();
        }

        if( aParam.has( "te012-mode" ) )
        {
            LPROG(lmclog,"Running with TE012 only.  Set parameter n-modes = 3");
            fTE012 = aParam["te012-mode"]().as_bool();
        }

        if( aParam.has( "te013-mode" ) )
        {
            LPROG(lmclog,"Running with TE013 only.  Set parameter n-modes = 4");
            fTE013 = aParam["te013-mode"]().as_bool();
        }


        if( aParam.has( "n-pixels" ) )
    	{
    		SetNPixels(aParam["n-pixels"]().as_int());
    	}

    	if( aParam.has( "plot-mode-maps" ) )
    	{
    		SetPlotModeMaps(aParam["plot-mode-maps"]().as_bool());
    	}

    	if ( aParam.has( "output-path" ) )
    	{
    		fOutputPath = aParam["output-path"]().as_string();
    	}

    	fAvgDotProductFactor.resize(fNModes);
    	for (unsigned m=0; m<fNModes; m++)
    	{
    		fAvgDotProductFactor[m].resize(fNModes);
        	for (unsigned n=0; n<fNModes; n++)
        	{
        		fAvgDotProductFactor[m][n].resize(fNModes);
        	}
    	}

    	return true;

    }

    int Field::GetNChannels()
    {
        return fNChannels;
    }
    void Field::SetNChannels( int aNumberOfChannels )
    {
     	fNChannels = aNumberOfChannels;
    }

    std::vector<std::vector<int>> Field::ModeSelect(bool bWaveguide, bool bNormCheck)
    {
    	int nModes = fNModes;
    	std::vector<std::vector<int>> tModeSet;
    	if ( !bNormCheck )
    	{
    	    if ( bWaveguide )
    	    {
    	    	tModeSet.push_back( {1,0,1,0} );
    	    }
    	    else
    	    {
    	    	if ( !fbMultiMode )
    	    	{
    	    	    if ( fTM111 )
    	    	    {
    	    	        tModeSet.push_back( {0,1,1,1} );
    	    	    }
    	    	    else if ( fTE012 )
    	    	    {
    	    	        tModeSet.push_back( {1,0,1,2} );
    	    	    }
    	    	    else if ( fTE013 )
    	    	    {
    	    	    	tModeSet.push_back( {1,0,1,3} );
    	    	    }
    	    	    else
    	    	    {
    	    	        tModeSet.push_back( {1,0,1,1} ); // default.
    	    	    }
    	    	}
    	    	else
    	    	{
    	    		tModeSet.push_back( {0,1,1,1} );
    	    		tModeSet.push_back( {1,0,1,1} );
    	    	}
    	    }
    	}
    	else
    	{
            for (int bTE=0; bTE<2; bTE++)
            {
                for (int l=0; l<nModes; l++)
                {
                    for (int m=1; m<nModes; m++)
                    {
                        for (int n=0; n<nModes; n++)
                        {
                        	tModeSet.push_back( {bTE,l,m,n} );
                        }
                    }
                }
            }
    	}

    	return tModeSet;
    }

    std::vector<std::vector<std::vector<std::vector<double>>>> Field::CalculateNormFactors(int nModes, bool bWaveguide)
    {

        LPROG(lmclog, "Calculating mode normalization factors ... " );

        std::vector<std::vector<std::vector<std::vector<double>>>> aModeNormFactor;
        aModeNormFactor.resize(2);

        for (int bTE=0; bTE<2; bTE++)
        {
            aModeNormFactor[bTE].resize(nModes);
            for (unsigned l=0; l<nModes; l++)
            {
                aModeNormFactor[bTE][l].resize(nModes);
                for (unsigned m=0; m<nModes; m++)
                {
                    aModeNormFactor[bTE][l][m].resize(nModes);
                    for (unsigned n=0; n<nModes; n++)
                    {
                        aModeNormFactor[bTE][l][m][n] = 0.;
                    }
                }
            }
        }

        std::vector<std::vector<int>> tModeSet = ModeSelect(bWaveguide, 0);

        for (int mu=0; mu<tModeSet.size(); mu++)
        {
            bool bTE = tModeSet[mu][0];
            int l = tModeSet[mu][1];
            int m = tModeSet[mu][2];
            int n = tModeSet[mu][3];

            aModeNormFactor[bTE][l][m][n] = 1./pow(Integrate(l,m,n,bTE,1),0.5);
        }

        return aModeNormFactor;
    }

    std::vector<std::vector<std::vector<std::vector<double>>>> Field::SetUnityNormFactors(int nModes, bool bWaveguide)
    {

        std::vector<std::vector<std::vector<std::vector<double>>>> aModeNormFactor;
        aModeNormFactor.resize(2);

        for (int bTE=0; bTE<2; bTE++)
        {
            aModeNormFactor[bTE].resize(nModes);
            for (unsigned l=0; l<nModes; l++)
            {
                aModeNormFactor[bTE][l].resize(nModes);
                for (unsigned m=0; m<nModes; m++)
                {
                    aModeNormFactor[bTE][l][m].resize(nModes);
                    for (unsigned n=0; n<nModes; n++)
                    {
                        aModeNormFactor[bTE][l][m][n] = 1.;
                    }
                }
            }
        }

        return aModeNormFactor;
    }


    void Field::CheckNormalization(int nModes, bool bWaveguide)
    {

        printf("\n \\int{|E_xlm|^2 dV} = \\mu / \\epsilon \\int{|H_xlm|^2 dV} ?\n\n");

        std::vector<std::vector<int>> tModeSet = ModeSelect(bWaveguide, 0);

        for (int mu=0; mu<tModeSet.size(); mu++)
        {
            bool bTE = tModeSet[mu][0];
            int l = tModeSet[mu][1];
            int m = tModeSet[mu][2];
            int n = tModeSet[mu][3];
            double normFactor = pow(GetNormFactors()[bTE][l][m][n],2.);

            if (bTE)
            {
            	if (!std::isnan(normFactor)&&(std::isfinite(normFactor)))
            	{
            	    printf("TE%d%d%d E %.4g H %.4g\n", l, m, n, Integrate(l,m,n,1,1)*normFactor,
            	    		LMCConst::MuNull()/LMCConst::EpsNull()*Integrate(l,m,n,1,0)*normFactor);
            	}
                else
                {
                    printf("TE%d%d%d is undefined.\n", l, m, n);
                }
            }
            else
            {
                if (!std::isnan(normFactor)&&(std::isfinite(normFactor)))
                {
                    printf("TM%d%d%d E %.4g H %.4g\n", l, m, n, Integrate(l,m,n,0,1)*normFactor,
                    		LMCConst::MuNull()/LMCConst::EpsNull()*Integrate(l,m,n,0,0)*normFactor);
                }
                else
                {
                    printf("TM%d%d%d is undefined.\n", l, m, n);
                }
            }
        }


        printf("\nThe modes normalized as above are available for use in the simulation.\n\n");
    }


    std::vector<std::vector<std::vector<std::vector<double>>>> Field::GetNormFactors()
    {
    	return fModeNormFactor;
    }

    void Field::SetNormFactors(std::vector<std::vector<std::vector<std::vector<double>>>> aNormFactor)
    {
    	fModeNormFactor = aNormFactor;
    }

    std::vector<std::vector<std::vector<double>>> Field::GetAvgDotProductFactor()
    {
    	return fAvgDotProductFactor;
    }

    void Field::SetAvgDotProductFactor(std::vector<std::vector<std::vector<double>>> aFactor)
    {
    	fAvgDotProductFactor = aFactor;
    }

    double Field::NormalizedEFieldMag(std::vector<double> field)
    {
    	double norm = 0;
    	auto it = field.begin();
    	while (it != field.end())
    	{
    		if (std::isfinite(*it)) norm += (*it)*(*it);
    		*it++;
    	}
    	return sqrt(norm);
    }

     bool FieldCore::Configure( const scarab::param_node& aParam )
     {
	return true;
     }

     double FieldCore::GetBesselNKZeros(int l, int m)
     {
     	return fBesselNKZeros[l][m];
     }

     double FieldCore::GetBesselNKPrimeZeros(int l, int m)
     {
     	return fBesselNKPrimeZeros[l][m];
     }

     void FieldCore::ReadBesselZeroes(std::string filename, bool prime)
     {
         std::ifstream input( filename );
         int n = 0; int k = 0; double zero = 0.;
         for( std::string line; getline( input, line ); )
         {
             std::stringstream ss(line);
             ss >> n;
             ss >> k;
             ss >> zero;
             if (prime)
             {
             	fBesselNKPrimeZeros.resize(n+1);
             	fBesselNKPrimeZeros[n].resize(k+1);
             	fBesselNKPrimeZeros[n][k] = zero;
             }
             else
             {
             	fBesselNKZeros.resize(n+1);
             	fBesselNKZeros[n].resize(k+1);
             	fBesselNKZeros[n][k] = zero;
             }

 //            printf("zero is %g and zero01 is %g\n\n", fInterface->fBesselNKPrimeZeros[n][k], fInterface->fBesselNKPrimeZeros[0][1]);
         }

     }


    double Field::GetCentralFrequency()
    {
    	return fCentralFrequency;
    }

    void Field::SetCentralFrequency( double aCentralFrequency )
    {
    	fCentralFrequency = aCentralFrequency;
    }

    int Field::GetNPixels()
    {
    	return fnPixels;
    }

    void Field::SetNPixels( int aNumberOfPixels )
    {
    	fnPixels = aNumberOfPixels;
    }

    int Field::GetNModes()
    {
    	return fNModes;
    }

    void Field::SetNModes( int aNumberOfModes )
    {
    	fNModes = aNumberOfModes;
    }

    double Field::GetDimX() const
    {
    	return fX;
    }

    void Field::SetDimX( double aDim )
    {
    	fX = aDim;
    }

    double Field::GetDimY() const
    {
    	return fY;
    }

    void Field::SetDimY( double aDim )
    {
    	fY = aDim;
    }

    double Field::GetDimR() const
    {
    	return fR;
    }

    void Field::SetDimR( double aDim )
    {
    	fR = aDim;
    }

    double Field::GetDimL() const
    {
    	return fL;
    }

    void Field::SetDimL( double aDim )
    {
    	fL = aDim;
    }

    double Field::GetCenterToShort() const
    {
    	return fCENTER_TO_SHORT;
    }

    void Field::SetCenterToShort( double aDistance )
    {
    	fCENTER_TO_SHORT = aDistance;
    }

    double Field::GetCenterToAntenna() const
    {
    	return fCENTER_TO_ANTENNA;
    }

    void Field::SetCenterToAntenna( double aDistance )
    {
    	fCENTER_TO_ANTENNA = aDistance;
    }

    bool Field::PlotModeMaps() const
    {
    	return fPlotModeMaps;
    }

    void Field::SetPlotModeMaps( bool aFlag )
    {
    	fPlotModeMaps = aFlag;
    }

    std::string Field::GetOutputPath()
    {
        return fOutputPath;
    }

    void Field::SetOutputPath( std::string aPath )
    {
     	fOutputPath = aPath;
    }


} /* namespace locust */

