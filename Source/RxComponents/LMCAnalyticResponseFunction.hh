/*
 * LMCAnalyticResponseFunction.hh
 *
 *  Created on: Jul 6, 2022
 *      Author: pslocum
 */

#ifndef LMCANALYTICRESPONSEFUNCTION_HH_
#define LMCANALYTICRESPONSEFUNCTION_HH_
#include "param.hh"
#include "LMCComplexFFT.hh"


namespace locust
{
 /*!
 @class LMCAnalyticResponseFunction
 @author P. Slocum
 @brief Base class to define analytic response functions
 @details
 Available configuration options:
 No input parameters
 */
    class AnalyticResponseFunction
    {

        public:
            AnalyticResponseFunction();
            virtual ~AnalyticResponseFunction();

            virtual bool Configure( const scarab::param_node& aNode );
            virtual bool GenerateTransferFunction() {return true;};
            virtual bool GenerateGreensFunction() {return true;};
            void SetGeneratingTF( bool aFlag );
            bool GetGeneratingTF();
            void SetInitialFreq( double aFreq );
            double GetInitialFreq();
            void SetTFarray( std::vector<std::complex<double>> aTFarray );
            std::vector<std::complex<double>> GetTFarray();

            void SetGFarray( std::vector <std::vector< std::vector< std::vector< std::vector<std::pair<double,std::pair<double,double> > > > > > > aGFarray );
            std::vector< std::vector<std::pair<double,std::pair<double,double> > > > GetGFarray( std::vector<std::vector<int>> );
            virtual void SetCavityQ( int bTE, int l, int m, int n, double aQ ) {};
            virtual double GetCavityQ(int bTE, int l, int m, int n ) {return 1.;};
            virtual void SetCavityFrequency(int bTE, int l, int m, int n,  double aFrequency ) {};
            virtual double GetCavityFrequency(int bTE, int l, int m, int n ) {return 0.;};
            virtual void SetDHOTimeResolution(int bTE, int l, int m, int n, double aTimeResolution ) {};
            virtual double GetDHOTimeResolution(int bTE, int l, int m, int n) {return 0.;};
            virtual void SetDHOThresholdFactor(int bTE, int l, int m, int n, double aThresholdFactor ) {};
            virtual double GetDHOThresholdFactor(int bTE, int l, int m, int n) {return 0.;};
            void SetNModes( int aNumberOfModes );
            int GetNModes();


        private:
            bool fGeneratingTF;
            std::vector<std::complex<double>> fTFarray;
            std::vector <std::vector< std::vector< std::vector< std::vector<std::pair<double,std::pair<double,double> > > > > > > fGFarray;
            double fInitialFreq;
            int fNModes;

};


} /* namespace locust */

#endif
