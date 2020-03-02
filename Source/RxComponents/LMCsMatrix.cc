/*
 * LMCSMatrix.cc
 *
 *  Created on: Feb 25, 2020
 *      Author: pslocum
 */

#include "LMCsMatrix.hh"
using std::string;


namespace locust
{

	LOGGER( lmclog, "SMatrix" );

    SMatrix::SMatrix():
    		fpatchImpedance( 100. ),
			fampImpedance( 50. ),
			fsMatrix2patch( 0 ),
			fsMatrix4patch( 0 ),
			fsMatrix6patch( 0 ),
			fsMatrix8patch( 0 )
    {
    }

    SMatrix::~SMatrix()
    {
    }

    bool SMatrix::Configure( const scarab::param_node& aParam )
    {

    	if( !PowerCombinerParent::Configure(aParam))
    	{
    		LERROR(lmclog,"Error configuring PowerCombiner class from SMatrix subclass");
    	}

    	SetVoltageDampingFactors();
    	return true;
    }

    std::vector<double> SMatrix::GetSmatrixElements()
    {
    	std::vector<double> sMatrix;

    	if (GetNElementsPerStrip() == 2) sMatrix = fsMatrix2patch;
    	if (GetNElementsPerStrip() == 4) sMatrix = fsMatrix4patch;
    	if (GetNElementsPerStrip() == 6) sMatrix = fsMatrix6patch;
    	if (GetNElementsPerStrip() == 8) sMatrix = fsMatrix8patch;

        return sMatrix;
    }


	bool SMatrix::SetVoltageDampingFactors()
	{
    	std::vector<double> sMatrix = GetSmatrixElements();
		for (unsigned z_index=0; z_index<GetNElementsPerStrip(); z_index++)
		{
			double aFactor = sMatrix[z_index+1]*sqrt(fampImpedance/fpatchImpedance);
			SetDampingFactor(z_index, aFactor);
		}

		return true;

	}




} /* namespace locust */
