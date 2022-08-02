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
            void SetGFarray( std::vector<std::pair<double,std::pair<double,double> > > aGFarray );
            std::vector<std::pair<double,std::pair<double,double> > > GetGFarray();


        private:
            bool fGeneratingTF;
            std::vector<std::complex<double>> fTFarray;
            std::vector<std::pair<double,std::pair<double,double> > > fGFarray;
            double fInitialFreq;








};


} /* namespace locust */

#endif
