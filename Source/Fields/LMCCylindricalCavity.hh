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
#include "LMCModeMapCavity.hh"
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
 @details  This class calculates field values induced by a moving electron in a
     cylindrical cavity.  The fields are calculated at either one probe location
     if the simulation parameter n-channels = 1, or at two probe locations if
     n-channels = 2.
 Available configuration options:
   - "cavity-radius" -- double [0.18] Radius of default cylindrical cavity, in meters.
   - "cavity-length" -- double [3.0] Length of default cylindrical cavity, in meters.
   - "cavity-probe-gain0" -- double [1.0] Gain of probe 0.
   - "cavity-probe-gain1" -- double [1.0] Gain of probe 1.
   - "cavity-probe-z0" -- double [0.0] Z location of probe 0, in meters.
   - "cavity-probe-z1" -- double [0.0] Z location of probe 1, in meters.
   - "cavity-probe-r-fraction0" -- double [0.5] Fractional radial position of probe 0 (unitless).
   - "cavity-probe-r-fraction1" -- double [0.5] Fractional radial position of probe 1 (unitless).
   - "cavity-probe-theta0" -- double [0.0] Azimuthal angle of probe 0 (radians).
   - "cavity-probe-theta1" -- double [0.0] Azimuthal angle of probe 1 (radians).
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
            virtual std::vector<double> GetNormalizedModeField(int l, int m, int n, std::vector<double> tKassParticleXP, bool includeOtherPols, bool teMode);
            virtual std::vector<double> GetTE_E(int l, int m, int n, double r, double theta, double z, bool includeOtherPols);
            virtual std::vector<double> GetTM_E(int l, int m, int n, double r, double theta, double z, bool includeOtherPols);
        	virtual double CalculateDotProductFactor(int l, int m, int n, std::vector<double> tKassParticleXP, std::vector<double> anE_normalized, double tThisEventNSamples);
            virtual double GetDotProductFactor(std::vector<double> tKassParticleXP, std::vector<double> anE_normalized, bool IntermediateFile);
            virtual bool InVolume(std::vector<double> tKassParticleXP);
            virtual void PrintModeMaps(int nModes, double zSlice, double thetaSlice);
            void PrintModeMapsLongSlice(int nModes, double thetaSlice);
            virtual std::vector<double> GetFieldAtProbe(int l, int m, int n, bool includeOtherPols, std::vector<double> tKassParticleXP, bool teMode);
            std::vector<double> GetCavityProbeZ();
            void SetCavityProbeZ ( double aZ, unsigned index );
            std::vector<double> GetCavityProbeR();
            void SetCavityProbeRFrac ( double aFraction, unsigned index );
            std::vector<double> GetCavityProbeGain();
            void SetCavityProbeTheta( double aTheta, unsigned index );
            std::vector<double> GetCavityProbeTheta();
            void SetCavityProbeGain( double aGain, unsigned index );


        private:
            bool fIntermediateFile;
            FieldCore* fFieldCore;
            std::vector<double> fCavityProbeZ;
            std::vector<double> fCavityProbeRFrac;
            std::vector<double> fCavityProbeTheta;
            std::vector<double> fProbeGain;
            bool fCaterpillarCavity;

    };



}; /* namespace locust */

#endif /* LMCCYLINDRICALCAVITY_HH_ */
