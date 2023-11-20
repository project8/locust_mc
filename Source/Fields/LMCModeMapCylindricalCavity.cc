/*
 * LMCModeMapCylindricalCavity.cc
 *
 *  Created on: Sep 29, 2023
 *      Author: pslocum
 */

#include "LMCModeMapCylindricalCavity.hh"


namespace locust
{

    LOGGER( lmclog, "ModeMapCylindricalCavity" );
    ModeMapCylindricalCavity::ModeMapCylindricalCavity():
        fModeMapTE_E(0.)
    {}

    ModeMapCylindricalCavity::~ModeMapCylindricalCavity(){}


    bool ModeMapCylindricalCavity::ReadModeMapTE_E(std::string aFilename)
    {
        std::fstream modeMapFile(aFilename.c_str(),std::ios::in);
        if (modeMapFile.fail())
        {
            LERROR(lmclog,"The mode map file \"" << aFilename <<"\" was not found.");
            return false;
        }

        while (!modeMapFile.eof())
        {
            std::string lineContent;
            while(std::getline(modeMapFile,lineContent))
            {
                std::vector<double> vectorLineContent;
                std::string token;
                std::stringstream ss(lineContent);
                int wordCount = 0;
                while (ss >> token)
                {
                	if (wordCount == 0) vectorLineContent.push_back(std::stod(token)); // radius position
                	else if (wordCount == 1) vectorLineContent.push_back(std::stod(token)); // theta position
                	else if (wordCount == 2) vectorLineContent.push_back(std::stod(token)); // mode E field value
                	else
                	{
                		LERROR(lmclog, "There are more columns than expected in the uploaded mode map file.");
                		return false;
                	}
                	++wordCount;
                }
                fModeMapTE_E.push_back(vectorLineContent);
                //printf("read r is %g, theta is %g, E is %g\n", fModeMapTE_E.back()[0], fModeMapTE_E.back()[1], fModeMapTE_E.back()[2]);
                //getchar();
            }
        }
        modeMapFile.close();

        return true;
    }

    std::vector<double> ModeMapCylindricalCavity::TE_E(double R, double L, int l, int m, int n, double r, double theta, double zKass, bool includeOtherPols)
    {
    	// This placeholder function does a coarse search through an uploaded mode map (stored in fModeMapTE_E)
    	// to find a data point near the requested location at (r,theta).  If it finds a suitable point, the function
    	// returns the value of the E field in the uploaded mode map at that point.  If it can't find
    	// a point, it returns 0.  This coarse grid search should be entirely replaced with something that returns
    	// an E-field value interpolated from the uploaded mode map.

    	// For testing, import the test mode map and inspect the output of this function like this:
    	//
    	// bin/LocustSim -c config/LocustCavityCCA.json "cavity-signal.upload-modemap-filename"="PozarExportTE011.txt" "cavity-signal.plot-mode-maps"=true
    	// root -l output/ModeMapOutput.root
    	// TE011_Etheta_z0->DrawCopy()
    	//
    	// The plotted mode map should appear as an exact copy of the uploaded mode map.
    	//
    	// Now, change the number of pixels used in the plotted mode map:
    	//
    	// bin/LocustSim -c config/LocustCavityCCA.json "cavity-signal.upload-modemap-filename"="PozarExportTE011.txt" "cavity-signal.plot-mode-maps"=true "cavity-signal.n-pixels"=201
    	// root -l output/ModeMapOutput.root
    	// TE011_Etheta_z0->DrawCopy()
    	//
    	// In the plot, we should see a few points that were extracted successfully from the uploaded mode map, but
    	// with most of the requested points missing.  This is due to the presently inadequate grid search.  With
    	// interpolation, we should be able to report reasonable E field values for any requested point over any number
    	// of pixels, and this should be seen as a smooth mode map in the output plot.


    	std::vector<double> TE_E;
    	if ((l==0)&&(m==1)&&(n==1)) // Hard-wire for lmn=011, which is the mode written to PozarExportTE011.txt .
    	{
    		for (int iTolerance=3; iTolerance>2; iTolerance--)
    		{
    			double tTolerance = pow(10.,-iTolerance);
    			for (unsigned j=0; j<fModeMapTE_E.size(); j++)
    			{
    				if ((fabs(1.-r/fModeMapTE_E[j][0]) < tTolerance) && (fabs(1.-theta/fModeMapTE_E[j][1]) < tTolerance))
    				{
    					// Found a near neighbor in the uploaded mode map:
    					printf("E at r=%g theta=%g is %g\n", r, theta, fModeMapTE_E[j][2]);
    					TE_E.push_back(fModeMapTE_E[j][2]);
    					return TE_E;  // Return the near neighbor.
    				}
    			}
    		}
    		TE_E.push_back(0.); // Never found a point.
    		return TE_E;  // Return 0.
    	}
    	else
    	{
    		TE_E.push_back(0.); // Wrong mode for right now.
    	}
        return TE_E;
    }

    std::vector<double> ModeMapCylindricalCavity::TE_H(double R, double L, int l, int m, int n, double r, double theta, double zKass, bool includeOtherPols)
    {
//    	LPROG( lmclog, "TE_H is presently not available as a mode map." );
    	std::vector<double> TE_H;
    	TE_H.push_back(0.);
    	return TE_H;
    }

    std::vector<double> ModeMapCylindricalCavity::TM_E(double R, double L, int l, int m, int n, double r, double theta, double zKass, bool includeOtherPols)
    {
//    	LPROG( lmclog, "TM_E is presently not available as a mode map." );
    	std::vector<double> TM_E;
    	TM_E.push_back(0.);
    	return TM_E;
    }

    std::vector<double> ModeMapCylindricalCavity::TM_H(double R, double L, int l, int m, int n, double r, double theta, double zKass, bool includeOtherPols)
    {
//    	LPROG( lmclog, "TM_H is presently not available as a mode map." );
    	std::vector<double> TM_H;
    	TM_H.push_back(0.);
        return TM_H;
    }

} /* namespace locust */

