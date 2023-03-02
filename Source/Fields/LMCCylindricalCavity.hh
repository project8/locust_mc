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
#include <boost/math/special_functions/bessel.hpp>
#include "LMCKassLocustInterface.hh"

#include <vector>

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


    class FieldCore
	{

    	public:

    	    FieldCore():
    	        fInterface( KLInterfaceBootstrapper::get_instance()->GetInterface() )
    	    {};
    	    virtual ~FieldCore(){};
            virtual std::vector<double> TE_E(double R, double L, int l, int m, int n, double r, double theta, double z, bool avgOverTheta){std::vector<double> x; return x;};
            virtual std::vector<double> TE_H(double R, double L, int l, int m, int n, double r, double theta, double z, bool avgOverTheta){std::vector<double> x; return x;};
            virtual std::vector<double> TM_E(double R, double L, int l, int m, int n, double r, double theta, double z, bool avgOverTheta){std::vector<double> x; return x;};
            virtual std::vector<double> TM_H(double R, double L, int l, int m, int n, double r, double theta, double z, bool avgOverTheta){std::vector<double> x; return x;};

    	private:
            kl_interface_ptr_t fInterface;


	};

    class PozarCylindrical: public FieldCore
    {
        public:
    	    PozarCylindrical():
    	        fInterface( KLInterfaceBootstrapper::get_instance()->GetInterface() )
            {};
    	    virtual ~PozarCylindrical(){};
            virtual std::vector<double> TE_E(double R, double L, int l, int m, int n, double r, double theta, double z, bool avgOverTheta);
            virtual std::vector<double> TE_H(double R, double L, int l, int m, int n, double r, double theta, double z, bool avgOverTheta);
            virtual std::vector<double> TM_E(double R, double L, int l, int m, int n, double r, double theta, double z, bool avgOverTheta);
            virtual std::vector<double> TM_H(double R, double L, int l, int m, int n, double r, double theta, double z, bool avgOverTheta);

        private:
            kl_interface_ptr_t fInterface;

    };



    class CylindricalCavity : public Field
    {

        public:
            CylindricalCavity();
            virtual ~CylindricalCavity();

            virtual bool Configure( const scarab::param_node& aParam);

            double Z_TE(int l, int m, int n, double fcyc) const;
            double Z_TM(int l, int m, int n, double fcyc) const;
            double Integrate(int l, int m, int n, bool teMode, bool eField);
            std::vector<double> GetDopplerFrequency(int l, int m, int n, std::vector<double> tKassParticleXP);
            std::vector<double> GetNormalizedModeField(int l, int m, int n, std::vector<double> tKassParticleXP);
            double GetDotProductFactor(std::vector<double> tKassParticleXP, std::vector<double> anE_normalized, bool IntermediateFile);

        private:
            kl_interface_ptr_t fInterface;
            FieldCore* fFieldCore;

    };



}; /* namespace locust */

#endif /* LMCCYLINDRICALCAVITY_HH_ */
