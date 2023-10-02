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
        	// Populate fModeMapTE_E with contents of file.
        	LPROG(lmclog,"Next step is to decide on a file format.");
        }
        modeMapFile.close();

    	return true;
    }

    std::vector<double> ModeMapCylindricalCavity::TE_E(double R, double L, int l, int m, int n, double r, double theta, double zKass, bool includeOtherPols)
    {
    	std::vector<double> TE_E;
    	// interpolate mode map modeMapTE_E at point r,theta,z, load it into TE_E.
        return TE_E;
    }

    std::vector<double> ModeMapCylindricalCavity::TE_H(double R, double L, int l, int m, int n, double r, double theta, double zKass, bool includeOtherPols)
    {
    	LPROG( lmclog, "TE_H is presently not available as a mode map." );
    	exit(-1);
    	std::vector<double> TE_H;
    	return TE_H;
    }

    std::vector<double> ModeMapCylindricalCavity::TM_E(double R, double L, int l, int m, int n, double r, double theta, double zKass, bool includeOtherPols)
    {
    	LPROG( lmclog, "TM_E is presently not available as a mode map." );
    	exit(-1);
    	std::vector<double> TM_E;
    	return TM_E;
    }

    std::vector<double> ModeMapCylindricalCavity::TM_H(double R, double L, int l, int m, int n, double r, double theta, double zKass, bool includeOtherPols)
    {
    	LPROG( lmclog, "TM_H is presently not available as a mode map." );
    	exit(-1);
    	std::vector<double> TM_H;
        return TM_H;
    }

} /* namespace locust */

