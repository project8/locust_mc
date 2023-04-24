/*
 * LMCCylindricalCavity.hh
 *
 *  Created on: Jun. 4, 2021
 *      Author: pslocum
 */

#ifndef LMCCYLINDRICALCAVITY_HH_
#define LMCCYLINDRICALCAVITY_HH_

#include "param.hh"

#include "logger.hh"
#include "LMCField.hh"
#include "LMCPozarCylindricalCavity.hh"
#include <boost/math/special_functions/bessel.hpp>

#include <vector>

#ifdef ROOT_FOUND
    #include "LMCRootHistoWriter.hh"
#endif


namespace locust
{
 /*!
 @class CylindricalCavity
 @author P. Slocum
 @brief Derived class to define CylindricalCavity fields.
 @details
 Available configuration options:
 No input parameters
 */


    class CylindricalCavity : public Field
    {

        public:
            CylindricalCavity();
            virtual ~CylindricalCavity();

            virtual bool Configure( const scarab::param_node& aParam);

            virtual double Z_TE(int l, int m, int n, double fcyc) const;
            virtual double Z_TM(int l, int m, int n, double fcyc) const;
            virtual double Integrate(int l, int m, int n, bool teMode, bool eField);
            virtual std::vector<double> GetDopplerFrequency(int l, int m, int n, std::vector<double> tKassParticleXP);
            virtual std::vector<double> GetNormalizedModeField(int l, int m, int n, std::vector<double> tKassParticleXP);
            virtual std::vector<std::vector<std::vector<double>>> CalculateNormFactors(int nModes, bool bTE);
            virtual std::vector<double> GetTE_E(int l, int m, int n, double r, double theta, double z, bool includeOtherPols);
            virtual double GetDotProductFactor(std::vector<double> tKassParticleXP, std::vector<double> anE_normalized, bool IntermediateFile);
            virtual void CheckNormalization(int nModes);
            virtual void PrintModeMaps(int nModes, bool bTE, double zSlice);
            virtual double GetFieldAtProbe(int l, int m, int n, bool includeOtherPols);
            double GetCavityProbeZ();
            void SetCavityProbeZ ( double aZ );
            double GetCavityProbeRFrac();
            void SetCavityProbeRFrac ( double aFraction );
            double GetCavityProbeGain();
            void SetCavityProbeGain( double aGain );


        private:
            FieldCore* fFieldCore;
            double fCavityProbeZ;
            double fCavityProbeRFrac;
            double fProbeGain;

    };



}; /* namespace locust */

#endif /* LMCCYLINDRICALCAVITY_HH_ */
