/*
 * LMCRectangularWaveguide.hh
 *
 *  Created on: Nov. 9, 2021
 *      Author: pslocum
 */

#ifndef LMCRECTANGULARWAVEGUIDE_HH_
#define LMCRECTANGULARWAVEGUIDE_HH_

#include "LMCField.hh"
#include "LMCKassLocustInterface.hh"


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


    class PozarRectangular: public FieldCore
    {
        public:
	        PozarRectangular():
	            fInterface( KLInterfaceBootstrapper::get_instance()->GetInterface() )
            {};
	        virtual ~PozarRectangular(){};
            virtual std::vector<double> TE_E(double dimX, double dimY, int m, int n, double xKass, double yKass, double fcyc) const;
            virtual std::vector<double> TE_H(double dimX, double dimY, int m, int n, double xKass, double yKass, double fcyc) const;
            virtual std::vector<double> TM_E(double dimX, double dimY, int m, int n, double xKass, double yKass, double fcyc) const;
            virtual std::vector<double> TM_H(double dimX, double dimY, int m, int n, double xKass, double yKass, double fcyc) const;

        private:
            kl_interface_ptr_t fInterface;

};



    class RectangularWaveguide : public Field
    {

        public:
            RectangularWaveguide();
            virtual ~RectangularWaveguide();

            virtual bool Configure( const scarab::param_node& aParam);

            double Z_TE(int l, int m, int n, double fcyc) const;
            double Z_TM(int l, int m, int n, double fcyc) const;
            double Integrate(int l, int m, int n, bool teMode, bool eField);
            std::vector<double> GetDopplerFrequency(int l, int m, int n, std::vector<double> tKassParticleXP);
            std::vector<double> GetNormalizedModeField(int l, int m, int n, std::vector<double> tKassParticleXP);
            double GetDotProductFactor(std::vector<double> tKassParticleXP, std::vector<double> aTE_E_normalized, bool IntermediateFile);




        private:
            double GetGroupVelocity(int m, int n, double fcyc);
            kl_interface_ptr_t fInterface;
            FieldCore* fFieldCore;


    };


}; /* namespace locust */

#endif /* LMCRECTANGULARWAVEGUIDE_HH_ */
