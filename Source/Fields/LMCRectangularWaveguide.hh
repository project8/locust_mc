/*
 * LMCRectangularWaveguide.hh
 *
 *  Created on: Nov. 9, 2021
 *      Author: pslocum
 */

#ifndef LMCRECTANGULARWAVEGUIDE_HH_
#define LMCRECTANGULARWAVEGUIDE_HH_

#include "LMCField.hh"
#include "LMCPozarRectangularWaveguide.hh"
#include "LMCKassLocustInterface.hh"
#ifdef ROOT_FOUND
    #include "LMCRootHistoWriter.hh"
#endif


#include <vector>

namespace locust
{
 /*!
 @class RectangularWaveguide
 @author P. Slocum
 @brief Derived class to define RectangularWaveguide fields.
 @details
 Available configuration options:
 No input parameters
 */


    class RectangularWaveguide : public Field
    {

        public:
            RectangularWaveguide();
            virtual ~RectangularWaveguide();

            virtual bool Configure( const scarab::param_node& aParam);

            virtual double Z_TE(int l, int m, int n, double fcyc) const;
            virtual double Z_TM(int l, int m, int n, double fcyc) const;
            virtual double Integrate(int l, int m, int n, bool teMode, bool eField);
            virtual std::vector<double> GetDopplerFrequency(int bTE, int l, int m, int n, std::vector<double> tKassParticleXP);
            virtual std::vector<double> GetNormalizedModeField(int l, int m, int n, std::vector<double> tKassParticleXP, bool includeOtherPols, bool bTE);
            virtual std::vector<std::vector<std::vector<double>>> CalculateNormFactors(int nModes, bool bTE);
            virtual double CalculateDotProductFactor(int bTE, int l, int m, int n, std::vector<double> tKassParticleXP, std::vector<double> anE_normalized, double tThisEventNSamples);
            virtual double GetDotProductFactor(std::vector<double> tKassParticleXP, std::vector<double> aTE_E_normalized, bool IntermediateFile);
            virtual void CheckNormalization(int nModes);
            virtual void PrintModeMaps(int nModes, bool bTE, double zSlice);
            double GetGroupVelocity(int m, int n, double fcyc);
            double ScaleEPoyntingVector(double fcyc);




        private:
            kl_interface_ptr_t fInterface;
            FieldCore* fFieldCore;


    };


}; /* namespace locust */

#endif /* LMCRECTANGULARWAVEGUIDE_HH_ */
