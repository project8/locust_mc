/*
 * LMCModeMapCavity.cc
 *
 *  Created on: Sep 29, 2023
 *      Author: pslocum
 */

#include "LMCModeMapCavity.hh"


namespace locust
{

    LOGGER( lmclog, "ModeMapCavity" );
    ModeMapCavity::ModeMapCavity():
        fModeMapTE_E(0.)
    {}

    ModeMapCavity::~ModeMapCavity(){}


    bool ModeMapCavity::ReadModeMapTE_E(std::string aFilename)
    {
        std::fstream modeMapFile(aFilename.c_str(),std::ios::in);
        if (modeMapFile.fail())
        {
            LERROR(lmclog,"The mode map file \"" << aFilename <<"\" was not found.");
            return false;
        }
        else
        {
            LWARN(lmclog,"Reading mode map file \"" << aFilename <<"\" .");
        }

	int nPixel1, nPixel2, nPixel3; //Represent number of pixels in provided Mode Map for each dimension. THESE NEED TO BE GENERALIZED!
	double init1 = 3.5e-5;
	double final1 = 0.006966;
	double init2 = 0.03142;
	double final2 = 6.253;
	double init3 = 0.016;
	double final3 = 0.085;
        nPixel1 = 100;
        nPixel2 = 100;
        nPixel3 = 3;
	fModeMapTE_E.resize(nPixel1);
        for(int i=0; i<nPixel1; i++)
	{
		fModeMapTE_E[i].resize(nPixel2);
		for(int j=0; j<nPixel2; j++)
		{
			fModeMapTE_E[i][j].resize(nPixel3);
		}
	}

        while (!modeMapFile.eof())
        {
            std::string lineContent;
            while(std::getline(modeMapFile,lineContent))
            {
		if(lineContent[0]=='#')
		{
			continue; //skips any comments
		}
                std::string token;
                std::stringstream ss(lineContent);
                int wordCount = 0;
		int i,j,k;
		double Etheta;
                while (ss >> token)
                {
                	if (wordCount == 0) i = (int)((std::stod(token)-init1)/(final1-init1)*nPixel1); // var1 position
                	else if (wordCount == 1) j = (int)((std::stod(token)-init2)/(final2-init2)*nPixel2); // var2 position
			else if (wordCount == 2) k = (int)((std::stod(token)-init3)/(final3-init3)*nPixel3); // var3 position
                	else if (wordCount == 3) Etheta = std::stod(token); // mode E field value
                	else
                	{
                		LERROR(lmclog, "There are more columns than expected in the uploaded mode map file.");
                		return false;
                	}
                	++wordCount;
                }
                fModeMapTE_E[i][j][k] = Etheta;
//		std::cout << "E read in at index " << i << " " << j << " " << k << ": " << fModeMapTE_E[i][j][k] << std::endl;
//                printf("read var1 is %g, var2 is %g, E is %g\n", fModeMapTE_E.back()[0], fModeMapTE_E.back()[1], fModeMapTE_E.back()[2]);
//                getchar();
            }
        }
        modeMapFile.close();

        return true;
    }

    std::vector<double> ModeMapCavity::TE_E(double dim1, double dim2, double dim3, int l, int m, int n, double var1, double var2, double zKass, bool includeOtherPols)
    {
    	// This placeholder function does a coarse search through an uploaded mode map (stored in fModeMapTE_E)
    	// to find a data point near the requested location at (var1,var2).  If it finds a suitable point, the function
    	// returns the value of the E field in the uploaded mode map at that point.  If it can't find
    	// a point, it returns 0.  This coarse grid search should be entirely replaced with something that returns
    	// an E-field value interpolated from the uploaded mode map.

    	// For testing, import the test mode map and inspect the output of this function like this:
    	//
    	// Cylindrical coordinates:
    	// bin/LocustSim -c config/LocustCavityCCA.json "cavity-signal.upload-modemap-filename"="PozarExport_Etheta_CylindricalTE011.txt" "cavity-signal.plot-mode-maps"=true
    	// root -l output/ModeMapOutput.root
    	// TE011_Etheta_z0mm->DrawCopy()
    	//
    	// Recangular coordinates:
    	// bin/LocustSim -c config/LocustCavityRectangular.json "cavity-signal.upload-modemap-filename"="PozarExport_Ex_RectangularTE011.txt" "cavity-signal.plot-mode-maps"=true
    	// root -l output/ModeMapOutput.root
    	// TE011_Ex_z0mm->DrawCopy()

    	// The plotted mode map should appear as an exact copy of the uploaded mode map.
    	//
    	// Now, change the number of pixels used in the plotted mode map:
    	//
    	// bin/LocustSim -c config/LocustCavityRectangular.json "cavity-signal.upload-modemap-filename"="PozarExport_Ex_RectangularTE011.txt" "cavity-signal.plot-mode-maps"=true "cavity-signal.n-pixels"=201
    	// root -l output/ModeMapOutput.root
    	// TE011_Ex_z0mm->DrawCopy()
    	//
    	// In the plot, we might see a few points that were extracted successfully from the uploaded mode map, but
    	// with most if not all of the requested points missing.  This is due to the presently inadequate grid search.  With
    	// interpolation, we should be able to report reasonable E field values for any requested point over any number
    	// of pixels, and this should be seen as a smooth and correct mode map in the output plot.

    	std::vector<double> TE_E;

    	// Below, we are presently hard-wired for lmn=011, which is the mode written to both available fake mode maps:
    	// A file containing a fake cylindrical TE011 mode map in r, theta, and E_theta:  PozarExport_Etheta_CylindricalTE011.txt,
    	// A file containing a fake rectangular TE011 mode map in x, y, and E_x:  PozarExport_Ex_RectangularTE011.txt .
    	//
    	if ((l==0)&&(m==1)&&(n==1))
    	{
		double var3 = zKass + 0.5*dim3;
		std::vector< int > CoordinateIndices = FindClosestCoordinate(var1, var2, var3, 0, dim1, fModeMapTE_E.size(), 0, dim2, fModeMapTE_E[0].size(), 0, dim3, fModeMapTE_E[0][0].size());
		std::vector< std::vector< int >> TetrahedronVertices = GetVerticesIndices(CoordinateIndices, var1, var2, var3, 0, dim1, fModeMapTE_E.size(), 0, dim2, fModeMapTE_E[0].size(), 0, dim3, fModeMapTE_E[0][0].size());
//		std::cout << "Indices: " << CoordinateIndices[0] << " " << CoordinateIndices[1] << " " << CoordinateIndices[2] << std::endl;
    		if(CoordinateIndices[0]!=0 or CoordinateIndices[1]!=0 or CoordinateIndices[2]!=0) 
    		{
    			// Found a near neighbor in the uploaded mode map:
    			//printf("E at var1=%g var2=%g var3=%g is %g\n", var1, var2, var3, fModeMapTE_E[CoordinateIndices[0]][CoordinateIndices[1]][CoordinateIndices[2]]);
			TE_E.push_back( InterpolateField(var1, var2, var3, TetrahedronVertices, 0, dim1, fModeMapTE_E.size(), 0, dim2, fModeMapTE_E[0].size(), 0, dim3, fModeMapTE_E[0][0].size() ));
    			return TE_E;  // Return the near neighbor.
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

    std::vector<double> ModeMapCavity::TE_H(double dim1, double dim2, double dim3, int l, int m, int n, double var1, double var2, double zKass, bool includeOtherPols)
    {
//    	LPROG( lmclog, "TE_H is presently not available as a mode map." );
    	std::vector<double> TE_H;
    	TE_H.push_back(0.);
    	return TE_H;
    }

    std::vector<double> ModeMapCavity::TM_E(double dim1, double dim2, double dim3, int l, int m, int n, double var1, double var2, double zKass, bool includeOtherPols)
    {
//    	LPROG( lmclog, "TM_E is presently not available as a mode map." );
    	std::vector<double> TM_E;
    	TM_E.push_back(0.);
    	return TM_E;
    }

    std::vector<double> ModeMapCavity::TM_H(double dim1, double dim2, double dim3, int l, int m, int n, double var1, double var2, double zKass, bool includeOtherPols)
    {
//    	LPROG( lmclog, "TM_H is presently not available as a mode map." );
    	std::vector<double> TM_H;
    	TM_H.push_back(0.);
        return TM_H;
    }


    std::vector< int > ModeMapCavity::FindClosestCoordinate(double var1, double var2, double var3, double dim1_min, double dim1_max, int dim1N, double dim2_min, double dim2_max, int dim2N, double dim3_min, double dim3_max, int dim3N)
    {
	//Finds coordinate indices with the floor of the index closest to that input variable for each dimension. Assumes a uniform grid in each of the 3 dimensions.
	std::vector< int > Coordinates(3);
	if( var1<dim1_min or var1>dim1_max or var2<dim2_min or var2>dim2_max or var3<dim3_min or var3>dim3_max )
	{
		Coordinates[0] = 0;
		Coordinates[1] = 0;
		Coordinates[2] = 0;
		return Coordinates;
	}
	else
	{
		Coordinates[0] = (int)((var1 - dim1_min)/(dim1_max - dim1_min)*dim1N + 0.5 - (var1<0)); //  "+ 0.5 - (var1<0)" means casting to int will round to the nearest int rather than truncate towards zero
		Coordinates[1] = (int)((var2 - dim2_min)/(dim2_max - dim2_min)*dim2N + 0.5 - (var1<0));
		Coordinates[2] = (int)((var3 - dim3_min)/(dim3_max - dim3_min)*dim3N + 0.5 - (var1<0));
		return Coordinates;
	}

    }

    std::vector< std::vector< int >> ModeMapCavity::GetVerticesIndices(std::vector<int> ClosestCoordinate, double var1, double var2, double var3, double dim1_min, double dim1_max, int dim1N, double dim2_min, double dim2_max, int dim2N, double dim3_min, double dim3_max, int dim3N)
    {
	//Assumes uniform grid spacing, but gives indicies for the coordinates of the 4 points from fModeMapTE_E that define the smallest tetrahedron enclosing a point in space
	double Closest_v1 = IndexToCoordinate(ClosestCoordinate[0], dim1_min, dim1_max, dim1N);
	double Closest_v2 = IndexToCoordinate(ClosestCoordinate[1], dim2_min, dim2_max, dim2N);
	double Closest_v3 = IndexToCoordinate(ClosestCoordinate[2], dim3_min, dim3_max, dim3N);
	std::vector< int > Vertex0 = ClosestCoordinate;
	std::vector< int > Vertex1 = ClosestCoordinate;
	std::vector< int > Vertex2 = ClosestCoordinate;
	std::vector< int > Vertex3 = ClosestCoordinate;
	if(var1>=Closest_v1)
	{
		Vertex1[0] += 1;
	}
	else
	{
		Vertex1[0] -= 1;
	}

        if(var2>=Closest_v2)
        {   
                Vertex2[1] += 1;
        }   
        else
        {   
                Vertex2[1] -= 1;
        } 

        if(var3>=Closest_v3)
        {   
                Vertex3[2] += 1;
        }   
        else
        {   
                Vertex3[2] -= 1;
        } 

	std::vector< std::vector< int >> VerticesIndices = {Vertex0, Vertex1, Vertex2, Vertex3};
	return VerticesIndices;

    }

    double ModeMapCavity::InterpolateField(double var1, double var2, double var3, std::vector< std::vector<int>> TetrahedronVertices, double dim1_min, double dim1_max, int dim1N, double dim2_min, double dim2_max, int dim2N, double dim3_min, double dim3_max, int dim3N)
    {
	//Does linear interpolation of a field at a point within a tetrahedron with known field at the vertices
	double x0, x1, x2, x3, y0, y1, y2, y3, z0, z1, z2, z3;
	x0 = IndexToCoordinate(TetrahedronVertices[0][0], dim1_min, dim1_max, dim1N);
        x1 = IndexToCoordinate(TetrahedronVertices[1][0], dim1_min, dim1_max, dim1N);
        x2 = IndexToCoordinate(TetrahedronVertices[2][0], dim1_min, dim1_max, dim1N);
        x3 = IndexToCoordinate(TetrahedronVertices[3][0], dim1_min, dim1_max, dim1N);
	y0 = IndexToCoordinate(TetrahedronVertices[0][1], dim2_min, dim2_max, dim2N);
        y1 = IndexToCoordinate(TetrahedronVertices[1][1], dim2_min, dim2_max, dim2N);
        y2 = IndexToCoordinate(TetrahedronVertices[2][1], dim2_min, dim2_max, dim2N);
        y3 = IndexToCoordinate(TetrahedronVertices[3][1], dim2_min, dim2_max, dim2N);
	z0 = IndexToCoordinate(TetrahedronVertices[0][2], dim3_min, dim3_max, dim3N);
        z1 = IndexToCoordinate(TetrahedronVertices[1][2], dim3_min, dim3_max, dim3N);
        z2 = IndexToCoordinate(TetrahedronVertices[2][2], dim3_min, dim3_max, dim3N);
        z3 = IndexToCoordinate(TetrahedronVertices[3][2], dim3_min, dim3_max, dim3N);
  for(int i=0; i<4; i++) //check if any indices are outside the grid of points used for interpolation
  {
	if(TetrahedronVertices[i][0]>=dim1N or TetrahedronVertices[i][2]>=dim3N) return 0.; //if r or z are outside of the grid, set field to zero
	if(TetrahedronVertices[i][0]<0 or TetrahedronVertices[i][2]<0) return 0.;

	if(TetrahedronVertices[i][1]>=dim2N) TetrahedronVertices[i][1] -= fModeMapTE_E[1].size(); //if theta is outside the the grid size wrap around back to 0 for periodicity
	if(TetrahedronVertices[i][1]<0) TetrahedronVertices[i][1] += fModeMapTE_E[1].size();
  }
  Eigen::MatrixXd m {
	{1., x0, y0, z0},
	{1., x1, y1, z1},
	{1., x2, y2, z2},
	{1., x3, y3, z3}
  };

  //std::cout << "m: " << m << std::endl; 

  Eigen::VectorXd v {{ fModeMapTE_E[TetrahedronVertices[0][0]][TetrahedronVertices[0][1]][TetrahedronVertices[0][2]], fModeMapTE_E[TetrahedronVertices[1][0]][TetrahedronVertices[1][1]][TetrahedronVertices[2][2]], fModeMapTE_E[TetrahedronVertices[2][0]][TetrahedronVertices[2][1]][TetrahedronVertices[2][2]], fModeMapTE_E[TetrahedronVertices[3][0]][TetrahedronVertices[3][1]][TetrahedronVertices[3][2]]}};

  //std::cout << "v: " << v << std::endl;

  Eigen::VectorXd Coef = m.inverse() * v;

  //std::cout << "Coef: " << Coef << std::endl;

  Eigen::VectorXd pos {{ 1., var1, var2, var3}};

  //std::cout << "pos: " << pos << std::endl;

  double interpolated_value = pos.dot(Coef);
  
  return interpolated_value;
  }

    double ModeMapCavity::IndexToCoordinate(int index, double min, double max, int nPixels)
    {
	return (double)index*(max - min)/((double)(nPixels-1)); //nPixels +1 since range in text file goes from bin midpoint to midpoint
    }

} /* namespace locust */

