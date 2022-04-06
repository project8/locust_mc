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


    class CylindricalCavity : public Field
    {

        public:
            CylindricalCavity();
            virtual ~CylindricalCavity();

            virtual bool Configure( const scarab::param_node& ) {return true;};

            std::vector<double> TE_E(int l, int m, int n, double r, double theta, double z, bool avgOverTheta) const;
            std::vector<double> TE_H(int l, int m, int n, double r, double theta, double z, bool avgOverTheta) const;
            std::vector<double> TM_E(int l, int m, int n, double r, double theta, double z, bool avgOverTheta) const;
            std::vector<double> TM_H(int l, int m, int n, double r, double theta, double z, bool avgOverTheta) const;
            double Z_TE(int l, int m, int n, double fcyc) const;
            double Z_TM(int l, int m, int n, double fcyc) const;
            double Integrate(int l, int m, int n, bool teMode, bool eField);
            double GetDopplerFrequency(int l, int m, int n, std::vector<double> tKassParticleXP);


        private:
            kl_interface_ptr_t fInterface;


    };


}; /* namespace locust */

#endif /* LMCCYLINDRICALCAVITY_HH_ */
