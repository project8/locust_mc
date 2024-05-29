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
 - "n-modes" : int -- [2] Range of l, m, and n indices used to configure available mode normalizations.
 Of the available normalized modes, modes to be included in the simulation itself are identified in the function
 LMCFieldCalculator::ModeSelect().
 - "n-pixels" : int -- [100] Number of pixels used in each dimension of the mode field definitions.
 - "plot-mode-maps": bool -- [false] Option to print all normalized (see n-modes above) mode maps to
 2D histograms in Root files, for inspection.
 */

    class FieldCore
    {

	    public:

    	    FieldCore();
    	    virtual ~FieldCore();

	    virtual bool Configure( const scarab::param_node& aParam );

	        // Cylindrical/rectangular Pozar cavities:
    	    // dim1 = r, dim2 = theta, dim3 = z
    	    // or
    	    // dim1 = x, dim2 = y, dim3 = z
            virtual std::vector<double> TE_E(double dim1, double dim2, double dim3, int l, int m, int n, double r, double theta, double z, bool includeOtherPols){return {0.};};
            virtual std::vector<double> TE_H(double dim1, double dim2, double dim3, int l, int m, int n, double r, double theta, double z, bool includeOtherPols){return {0.};};
            virtual std::vector<double> TM_E(double dim1, double dim2, double dim3, int l, int m, int n, double r, double theta, double z, bool includeOtherPols){return {0.};};
            virtual std::vector<double> TM_H(double dim1, double dim2, double dim3, int l, int m, int n, double r, double theta, double z, bool includeOtherPols){return {0.};};

            // Rectangular Pozar waveguide:
            virtual std::vector<double> TE_E(double dimX, double dimY, int m, int n, double xKass, double yKass, double fcyc){return {0.};};
            virtual std::vector<double> TE_H(double dimX, double dimY, int m, int n, double xKass, double yKass, double fcyc){return {0.};};
            virtual std::vector<double> TM_E(double dimX, double dimY, int m, int n, double xKass, double yKass, double fcyc){return {0.};};
            virtual std::vector<double> TM_H(double dimX, double dimY, int m, int n, double xKass, double yKass, double fcyc){return {0.};};

            // Imported mode map:
            virtual bool ReadModeMapTE_E(std::string aFilename){return 0;};
            virtual bool ReadModeMapTE_H(std::string aFilename){return 0;};
            virtual bool ReadModeMapTM_E(std::string aFilename){return 0;};
            virtual bool ReadModeMapTM_H(std::string aFilename){return 0;};

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
            virtual std::vector<double> GetNormalizedModeField(int l, int m, int n, std::vector<double> tKassParticleXP, bool includeOtherPols, bool teMode) {return {0.};};
            double NormalizedEFieldMag(std::vector<double> field);
            virtual std::vector<std::vector<std::vector<double>>> CalculateNormFactors(int nModes, bool bTE) {return {{{0.}}};};
            std::vector<std::vector<std::vector<double>>> SetUnityNormFactors(int nModes);
            virtual std::vector<double> GetTE_E(int l, int m, int n, double r, double theta, double z, bool includeOtherPols) {return {0.};};
            virtual std::vector<double> GetTM_E(int l, int m, int n, double r, double theta, double z, bool includeOtherPols) {return {0.};};
            virtual double CalculateDotProductFactor(int l, int m, int n, std::vector<double> tKassParticleXP, std::vector<double> aTE_E_normalized, double tThisEventNSamples) {return {0.};};
            virtual double GetDotProductFactor(std::vector<double> tKassParticleXP, std::vector<double> aTE_E_normalized, bool IntermediateFile) {return {0.};};
            virtual bool InVolume(std::vector<double> tKassParticleXP){return false;};
            virtual void CheckNormalization(int nModes){};
            virtual void PrintModeMaps(int nModes, double zSlice, double thetaSlice){};
            virtual std::vector<double> GetFieldAtProbe(int l, int m, int n, bool includeOtherPols, std::vector<double> tKassParticleXP, bool teMode){return {0.};};
            virtual double ScaleEPoyntingVector(double fcyc){return 0.;};

            std::vector<std::vector<std::vector<double>>> GetNormFactorsTE();
            void SetNormFactorsTE(std::vector<std::vector<std::vector<double>>> aNormFactor);
            std::vector<std::vector<std::vector<double>>> GetNormFactorsTM();
            void SetNormFactorsTM(std::vector<std::vector<std::vector<double>>> aNormFactor);
            std::vector<std::vector<std::vector<double>>> GetAvgDotProductFactor();
            void SetAvgDotProductFactor(std::vector<std::vector<std::vector<double>>> aFactor);
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
            double GetCenterToShort() const;
            void SetCenterToShort( double aDistance );
            double GetCenterToAntenna() const;
            void SetCenterToAntenna( double aDistance );
            bool PlotModeMaps() const;
            void SetPlotModeMaps( bool aFlag );
            void SetOutputPath( std::string aPath );
            std::string GetOutputPath();



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
            double fCENTER_TO_SHORT;
            double fCENTER_TO_ANTENNA;
            std::vector<std::vector<std::vector<double>>> fAvgDotProductFactor;
            bool fPlotModeMaps;
            std::string fOutputPath;

    };


}; /* namespace locust */

#endif /* LMCFIELD_HH_ */
