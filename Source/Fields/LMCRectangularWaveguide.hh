/*
 * LMCRectangularWaveguide.hh
 *
 *  Created on: Nov. 9, 2021
 *      Author: pslocum
 */

#ifndef LMCRECTANGULARWAVEGUIDE_HH_
#define LMCRECTANGULARWAVEGUIDE_HH_

#include "param.hh"

#include "logger.hh"
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


    class RectangularWaveguide : public Field
    {

        public:
            RectangularWaveguide();
            virtual ~RectangularWaveguide();

            virtual bool Configure( const scarab::param_node& ){return true;};

            std::vector<double> TE_E(int m, int n, double xKass, double yKass, double fcyc) const;
            std::vector<double> TE_H(int m, int n, double xKass, double yKass, double fcyc) const;
            std::vector<double> TM_E(int m, int n, double xKass, double yKass, double fcyc) const;
            std::vector<double> TM_H(int m, int n, double xKass, double yKass, double fcyc) const;
            double Z_TE(int l, int m, int n, double fcyc) const;
            double Z_TM(int l, int m, int n, double fcyc) const;
            double Integrate(int l, int m, int n, bool teMode, bool eField);
            std::vector<double> GetDopplerFrequency(int l, int m, int n, std::vector<double> tKassParticleXP);


        private:
            double GetGroupVelocity(int m, int n, double fcyc);
            kl_interface_ptr_t fInterface;


    };


}; /* namespace locust */

#endif /* LMCRECTANGULARWAVEGUIDE_HH_ */
