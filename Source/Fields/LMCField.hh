/*
 * LMCField.hh
 *
 *  Created on: Jun. 4, 2021
 *      Author: pslocum
 */

#ifndef LMCFIELD_HH_
#define LMCFIELD_HH_

#include "param.hh"
#include "logger.hh"
#include "LMCConst.hh"


#include <vector>

namespace locust
{
 /*!
 @class Field
 @author P. Slocum
 @brief Base class to characterize Field selection
 @details
 Available configuration options:
 No input parameters
 */

    class FieldCore
    {

	    public:

    	    FieldCore();
    	    virtual ~FieldCore();

	        // Cylindrical cavity:
            virtual std::vector<double> TE_E(double R, double L, int l, int m, int n, double r, double theta, double z, bool avgOverTheta){return {0.};};
            virtual std::vector<double> TE_H(double R, double L, int l, int m, int n, double r, double theta, double z, bool avgOverTheta){return {0.};};
            virtual std::vector<double> TM_E(double R, double L, int l, int m, int n, double r, double theta, double z, bool avgOverTheta){return {0.};};
            virtual std::vector<double> TM_H(double R, double L, int l, int m, int n, double r, double theta, double z, bool avgOverTheta){return {0.};};

            // Rectangular waveguide:
            virtual std::vector<double> TE_E(double dimX, double dimY, int m, int n, double xKass, double yKass, double fcyc){return {0.};};
            virtual std::vector<double> TE_H(double dimX, double dimY, int m, int n, double xKass, double yKass, double fcyc){return {0.};};
            virtual std::vector<double> TM_E(double dimX, double dimY, int m, int n, double xKass, double yKass, double fcyc){return {0.};};
            virtual std::vector<double> TM_H(double dimX, double dimY, int m, int n, double xKass, double yKass, double fcyc){return {0.};};

            void ReadBesselZeroes(std::string filename, bool prime);
            double GetBesselNKZeros(int l, int m);
            double GetBesselNKPrimeZeros(int l, int m);

        private:
            std::vector<std::vector<double> > fBesselNKZeros, fBesselNKPrimeZeros;


    };


    class Field
    {

        public:
            Field();
            virtual ~Field();

            virtual bool Configure( const scarab::param_node& aParam );

            virtual double Z_TM(int l, int m, int n, double fcyc) const {return {0.};};
            virtual double Z_TE(int l, int m, int n, double fcyc) const {return {0.};};
            virtual double Integrate(int l, int m, int n, bool teMode, bool eField){return 0.;};
            virtual std::vector<double> GetDopplerFrequency(int l, int m, int n, std::vector<double> tKassParticleXP) {return {0.};};
            virtual std::vector<double> GetNormalizedModeField(int l, int m, int n, std::vector<double> tKassParticleXP) {return {0.};};
            virtual std::vector<std::vector<std::vector<double>>> CalculateNormFactors(int nModes, bool bTE) {return {{{0.}}};};
            virtual std::vector<double> GetTE_E(int l, int m, int n, double r, double theta, double z, bool avgOverTheta) {return {0.};};
            virtual double GetDotProductFactor(std::vector<double> tKassParticleXP, std::vector<double> aTE_E_normalized, bool IntermediateFile) {return {0.};};
            virtual void CheckNormalization(int nModes){};
            virtual void PrintModeMaps(int nModes, bool bTE, double zSlice){};
            virtual double GetFieldAtProbe(int l, int m, int n, bool avgOverTheta){return 0.;};
            virtual double ScaleEPoyntingVector(double fcyc){return 0.;};


            std::vector<std::vector<std::vector<double>>> GetNormFactorsTE();
            void SetNormFactorsTE(std::vector<std::vector<std::vector<double>>> aNormFactor);
            std::vector<std::vector<std::vector<double>>> GetNormFactorsTM();
            void SetNormFactorsTM(std::vector<std::vector<std::vector<double>>> aNormFactor);
            double GetCentralFrequency();
            void SetCentralFrequency( double aCentralFrequency );
            int GetNPixels();
            void SetNPixels( int aNumberOfPixels );
            int GetNModes();
            void SetNModes( int aNumberOfModes );
            double GetDimX() const;
            void SetDimX( double aDim );
            double GetDimY() const;
            void SetDimY( double aDim );
            double GetDimR() const;
            void SetDimR( double aDim );
            double GetDimL() const;
            void SetDimL( double aDim );



        private:
            int fNModes;
            std::vector<std::vector<std::vector<double>>> fModeNormFactorTE;  // 3D vector [n-modes][n-modes][n-modes].
            std::vector<std::vector<std::vector<double>>> fModeNormFactorTM;  // 3D vector [n-modes][n-modes][n-modes].
            double fCentralFrequency;
            int fnPixels;
            double fR;  // Cylindrical cavity dimenions.
            double fL;
            double fX;  // Rectangular waveguide dimensions.
            double fY;




    };


}; /* namespace locust */

#endif /* LMCFIELD_HH_ */
