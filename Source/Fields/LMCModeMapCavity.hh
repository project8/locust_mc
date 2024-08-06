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
#include <Eigen/Dense>

namespace locust
{
 /*!
 @class ModeMapCavity
 @author P. Slocum
 @brief Derived class to define Cavity fields from an uploaded mode map file.
 @details
 Usage:
     	// Cylindrical coordinates:
    	// bin/LocustSim -c config/LocustCavityCCA_ModeMapTest.json "cavity-signal.upload-modemap-filename"="fieldsExportTest.fld" "cavity-signal.plot-mode-maps"=true
    	// root -l output/ModemapOutput.root
    	// TE011_Etheta_z0mm->DrawCopy()
    	//
    	// Recangular coordinates:
    	// bin/LocustSim -c config/LocustCavityRectangular.json "cavity-signal.upload-modemap-filename"="[putFilenameHere]" "cavity-signal.plot-mode-maps"=true
    	// root -l output/ModeMapOutput.root
    	// TE011_Ex_z0mm->DrawCopy()

    	// The plotted mode map should appear as an exact copy of the uploaded mode map.

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
            virtual bool ReadModeMapTE_E(std::string aFilename, const scarab::param_node& aParam);
            std::vector< int > FindClosestCoordinate(double var1, double var2, double var3);
            std::vector< std::vector< int >> GetVerticesIndices(std::vector<int> ClosestCoordinate, double var1, double var2, double var3);
            double InterpolateField(double var1, double var2, double var3, std::vector< std::vector<int>> TetrahedronVertices, int component);
            double IndexToCoordinate(int index, double min, double max, int nPixels);
            void MatchCavityDimensions(const scarab::param_node& aParam);

        private:
            int fnPixel1, fnPixel2, fnPixel3;
            double fDim1_min, fDim2_min, fDim3_min;
            double fDim1_max, fDim2_max, fDim3_max;
            std::vector< std::vector< std::vector< std::vector< double >>>> fModeMapTE_E;

    };


}; /* namespace locust */

#endif /* LMCMODEMAPCAVITY_HH_ */
