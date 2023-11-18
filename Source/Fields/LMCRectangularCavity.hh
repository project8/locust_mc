/*
 * LMCRectangularCavity.hh
 *
 *  Created on: Nov. 16, 2023
 *      Author: pslocum
 */

#ifndef LMCRECTANGULARCAVITY_HH_
#define LMCRECTANGULARCAVITY_HH_

#include "param.hh"

#include "logger.hh"
#include "LMCField.hh"
#include "LMCPozarRectangularCavity.hh"
//#include "LMCModeMapRectangularCavity.hh"

#include <vector>

#ifdef ROOT_FOUND
    #include "LMCRootHistoWriter.hh"
#endif


namespace locust
{
 /*!
 @class RectangularCavity
 @author P. Slocum
 @brief Derived class to define RectangularCavity fields.
 @details  This class calculates field values induced by a moving electron in a
     rectangular cavity.  The fields are calculated at either one probe location
     if the simulation parameter n-channels = 1, or at two probe locations if
     n-channels = 2.
 Available configuration options:
 */


    class RectangularCavity : public Field
    {

        public:
            RectangularCavity();
            virtual ~RectangularCavity();

            virtual bool Configure( const scarab::param_node& aParam);

            virtual double Z_TE(int l, int m, int n, double fcyc) const;
            virtual double Z_TM(int l, int m, int n, double fcyc) const;
            virtual double Integrate(int l, int m, int n, bool teMode, bool eField);
            virtual std::vector<double> GetDopplerFrequency(int l, int m, int n, std::vector<double> tKassParticleXP);
            virtual std::vector<double> GetNormalizedModeField(int l, int m, int n, std::vector<double> tKassParticleXP, bool includeOtherPols, bool teMode);
            virtual std::vector<std::vector<std::vector<double>>> CalculateNormFactors(int nModes, bool bTE);
            virtual std::vector<double> GetTE_E(int l, int m, int n, double r, double theta, double z, bool includeOtherPols);
            virtual std::vector<double> GetTM_E(int l, int m, int n, double r, double theta, double z, bool includeOtherPols);
        	virtual double CalculateDotProductFactor(int l, int m, int n, std::vector<double> tKassParticleXP, std::vector<double> anE_normalized, double tThisEventNSamples);
            virtual double GetDotProductFactor(std::vector<double> tKassParticleXP, std::vector<double> anE_normalized, bool IntermediateFile);
            virtual bool InVolume(std::vector<double> tKassParticleXP);
            virtual void CheckNormalization(int nModes);
            virtual void PrintModeMaps(int nModes, double zSlice);
            void PrintModeMapsLongSlice(int nModes, double xSlice);
            virtual std::vector<double> GetFieldAtProbe(int l, int m, int n, bool includeOtherPols, std::vector<double> tKassParticleXP, bool teMode);
            std::vector<double> GetCavityProbeZ();
            void SetCavityProbeZ ( double aZ, unsigned index );
            std::vector<double> GetCavityProbeRFrac();
            void SetCavityProbeRFrac ( double aFraction, unsigned index );
            std::vector<double> GetCavityProbeGain();
            void SetCavityProbeTheta( double aTheta, unsigned index );
            std::vector<double> GetCavityProbeTheta();
            void SetCavityProbeGain( double aGain, unsigned index );

        private:
            FieldCore* fFieldCore;
            std::vector<double> fCavityProbeZ;
            std::vector<double> fCavityProbeRFrac;
            std::vector<double> fCavityProbeTheta;
            std::vector<double> fProbeGain;

    };



}; /* namespace locust */

#endif /* LMCRECTANGULARCAVITY_HH_ */
