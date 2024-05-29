/*
 * LMCModeMapCavity.hh
 *
 *  Created on: Sept. 29, 2023
 *      Author: pslocum
 */

#ifndef LMCMODEMAPCAVITY_HH_
#define LMCMODEMAPCAVITY_HH_

#include "param.hh"

#include "logger.hh"
#include "LMCField.hh"

#include <vector>
#include <Dense>

namespace locust
{
 /*!
 @class ModeMapCavity
 @author P. Slocum
 @brief Derived class to define Cavity fields from an uploaded mode map file.
 @details
 Available configuration options:
 No input parameters
 */


    class ModeMapCavity: public FieldCore
    {
        public:

    	    ModeMapCavity();
		    virtual ~ModeMapCavity();
            virtual std::vector<double> TE_E(double dim1, double dim2, double dim3, int l, int m, int n, double r, double theta, double z, bool includeOtherPols);
            virtual std::vector<double> TE_H(double dim1, double dim2, double dim3, int l, int m, int n, double r, double theta, double z, bool includeOtherPols);
            virtual std::vector<double> TM_E(double dim1, double dim2, double dim3, int l, int m, int n, double r, double theta, double z, bool includeOtherPols);
            virtual std::vector<double> TM_H(double dim1, double dim2, double dim3, int l, int m, int n, double r, double theta, double z, bool includeOtherPols);

	    virtual bool Configure(const scarab::param_node& aParam);
            virtual bool ReadModeMapTE_E(std::string aFilename);
	    std::vector< int > FindClosestCoordinate(double var1, double var2, double var3);
	    std::vector< std::vector< int >> GetVerticesIndices(std::vector<int> ClosestCoordinate, double var1, double var2, double var3); 
	    double InterpolateField(double var1, double var2, double var3, std::vector< std::vector<int>> TetrahedronVertices, int component);
	    double IndexToCoordinate(int index, double min, double max, int nPixels);
        private:
	    int nPixel1, nPixel2, nPixel3;
	    double dim1_min, dim2_min, dim3_min;
	    double dim1_max, dim2_max, dim3_max;
            std::vector< std::vector< std::vector< std::vector< double >>>> fModeMapTE_E;

    };


}; /* namespace locust */

#endif /* LMCMODEMAPCAVITY_HH_ */
