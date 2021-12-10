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


    class Field
    {

        public:
            Field();
            virtual ~Field();

            bool Configure( const scarab::param_node& aNode );


            // size of field vectors will be number of components in field value at (r,theta,z)

            // cylindrical cavity
            virtual std::vector<double> TE_E(int l, int m, int n, double r, double theta, double z, double fcyc) const {return {0.};};
            virtual std::vector<double> TE_H(int l, int m, int n, double r, double theta, double z, double fcyc) const {return {0.};};
            virtual std::vector<double> TM_E(int l, int m, int n, double r, double theta, double z, double fcyc) const {return {0.};};
            virtual std::vector<double> TM_H(int l, int m, int n, double r, double theta, double z, double fcyc) const {return {0.};};
            virtual double Z_TM(int l, int m, int n) const {return {0.};};
            virtual double Z_TE(int l, int m, int n) const {return {0.};};


            // rectangular waveguide
            virtual std::vector<double> TE_E(int m, int n, double x, double y, double fcyc) const {return {0.};};
            virtual std::vector<double> TE_H(int m, int n, double x, double y, double fcyc) const {return {0.};};
            virtual std::vector<double> TM_E(int m, int n, double x, double y, double fcyc) const {return {0.};};
            virtual std::vector<double> TM_H(int m, int n, double x, double y, double fcyc) const {return {0.};};


            virtual double Integrate(int l, int m, int n, bool teMode, bool eField){return 0.;};

            std::vector<std::vector<std::vector<double>>> GetNormFactorsTE();
            void SetNormFactorsTE(std::vector<std::vector<std::vector<double>>> aNormFactor);
            std::vector<std::vector<std::vector<double>>> GetNormFactorsTM();
            void SetNormFactorsTM(std::vector<std::vector<std::vector<double>>> aNormFactor);
            double GetCentralFrequency();
            void SetCentralFrequency( double aCentralFrequency );

        private:
            std::vector<std::vector<std::vector<double>>> fModeNormFactorTE;  // 3D vector [n-modes][n-modes][n-modes].
            std::vector<std::vector<std::vector<double>>> fModeNormFactorTM;  // 3D vector [n-modes][n-modes][n-modes].
            double fCentralFrequency;

    };


}; /* namespace locust */

#endif /* LMCFIELD_HH_ */
