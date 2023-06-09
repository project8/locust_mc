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
        fR( 0.18 ),
        fL( 3.0 ),
        fX( 0.010668 ),
        fY( 0.004318 )
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

    	if( aParam.has( "central-frequency" ) )
    	{
    		fCentralFrequency= 2.*LMCConst::Pi()*aParam["central-frequency"]().as_double();
    	}
    	if( aParam.has( "n-pixels" ) )
    	{
    		SetNPixels(aParam["n-pixels"]().as_int());
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


    std::vector<std::vector<std::vector<double>>> Field::GetNormFactorsTE()
    {
    	return fModeNormFactorTE;
    }

    void Field::SetNormFactorsTE(std::vector<std::vector<std::vector<double>>> aNormFactor)
    {
    	fModeNormFactorTE = aNormFactor;
    }

    std::vector<std::vector<std::vector<double>>> Field::GetNormFactorsTM()
    {
    	return fModeNormFactorTM;
    }

    void Field::SetNormFactorsTM(std::vector<std::vector<std::vector<double>>> aNormFactor)
    {
    	fModeNormFactorTM = aNormFactor;
    }

    std::vector<std::vector<std::vector<double>>> Field::GetAvgDotProductFactor()
    {
    	return fAvgDotProductFactor;
    }

    void Field::SetAvgDotProductFactor(std::vector<std::vector<std::vector<double>>> aFactor)
    {
    	fAvgDotProductFactor = aFactor;
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





} /* namespace locust */

