/*
 * LMCSMatrix.hh
 *
 *  Created on: Feb 28, 2020
 *      Author: pslocum
 */

#ifndef LMCSMATRIX_HH_
#define LMCSMATRIX_HH_

#include "LMCPowerCombinerParent.hh"
#include "param.hh"
#include "logger.hh"
#include <iostream>

namespace locust
{
 /*!
 @class SMatrix
 @author P. Slocum
 @brief Derived class describing a power combiner in terms of its S-matrix from HFSS.
 @details
 Available configuration options:
 No input parameters
 */


    class SMatrix: public PowerCombinerParent
    {

        public:
            SMatrix();
            virtual ~SMatrix();
            virtual bool Configure( const scarab::param_node& aNode );
        	virtual bool SetVoltageDampingFactors();



        private:

        	std::vector<double> GetSmatrixElements();
        	double fpatchImpedance;
        	double fampImpedance;

        	// 7/8 power combiner S-matrices (traveling wave configuration) from HFSS:
            std::vector<double> fsMatrix2patch = {0.12, 0.43, 0.43};
            std::vector<double> fsMatrix4patch = {0.12, 0.3, 0.43, 0.43, 0.3};
            std::vector<double> fsMatrix6patch = {0.12, 0.24, 0.3, 0.43, 0.43, 0.3, 0.24};
            std::vector<double> fsMatrix8patch = {0.12, 0.17, 0.24, 0.3, 0.43, 0.43, 0.3, 0.24, 0.17};
            // end 7/8 combiner S-matrices.




    };


} /* namespace locust */

#endif /* LMCSMATRIX_HH_ */
