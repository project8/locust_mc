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
        fModeMapTE_E(0.),
        fDim1_min(0.),
        fDim1_max(0.0762),
        fDim2_min(0.),
        fDim2_max(6.284),
        fDim3_min(0.),
        fDim3_max(0.1524),
        fnPixel1(10),
        fnPixel2(10),
        fnPixel3(10)
    {}

    ModeMapCavity::~ModeMapCavity(){}

    bool ModeMapCavity::Configure( const scarab::param_node& aParam)
    {
        if( aParam.has( "n-pixel1" ) )
        {
            fnPixel1 =  aParam["n-pixel1"]().as_int();
        }

	if( aParam.has( "dim1-min" ) )
        {
            fDim1_min =  aParam["dim1-min"]().as_double();
        }

        if( aParam.has( "dim1-max" ) )
        {
            fDim1_max =  aParam["dim1-max"]().as_double();
        }

        if( aParam.has( "n-pixel2" ) )
        {   
            fnPixel2 =  aParam["n-pixel2"]().as_int();
        }   

        if( aParam.has( "dim2-min" ) )
        {
            fDim2_min =  aParam["dim2-min"]().as_double();
        }

        if( aParam.has( "dim2-max" ) )
        {
            fDim2_max =  aParam["dim2-max"]().as_double();
        }

        if( aParam.has( "n-pixel3" ) )
        {   
            fnPixel3 =  aParam["n-pixel3"]().as_int();
        }   

        if( aParam.has( "dim3-min" ) )
        {
            fDim3_min =  aParam["dim3-min"]().as_double();
        }

        if( aParam.has( "dim3-max" ) )
        {
            fDim3_max =  aParam["dim3-max"]().as_double();
        }

	return true;
    }

    bool ModeMapCavity::ReadModeMapTE_E(std::string aFilename, const scarab::param_node& aParam)
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

        fModeMapTE_E.resize(fnPixel1);
        for(int i=0; i<fnPixel1; i++)
        {
            fModeMapTE_E[i].resize(fnPixel2);
            for(int j=0; j<fnPixel2; j++)
            {
                fModeMapTE_E[i][j].resize(fnPixel3);
                for(int k=0; k<fnPixel3; k++)
                {
                    fModeMapTE_E[i][j][k].resize(3);
                }
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
                double Erho,Etheta,Ez;
                while (ss >> token)
                {
                    if (wordCount == 0) i = (int)((std::stod(token)-fDim1_min)/(fDim1_max-fDim1_min)*(fnPixel1)); // var1 position
                    else if (wordCount == 1) j = (int)((std::stod(token)-fDim2_min)/(fDim2_max-fDim2_min)*(fnPixel2)); // var2 position
                    else if (wordCount == 2) k = (int)((std::stod(token)-fDim3_min)/(fDim3_max-fDim3_min)*(fnPixel3)); // var3 position
                    else if (wordCount == 3) Erho = std::stod(token); // mode E field value
                    else if (wordCount == 4) Etheta = std::stod(token); // mode E field value
                    else if (wordCount == 5) Ez = std::stod(token); // mode E field value
                    else
                    {
                        LERROR(lmclog, "There are more columns than expected in the uploaded mode map file.");
                        return false;
                    }
                    ++wordCount;
                }

                if ((i>=fnPixel1) or (j>=fnPixel2) or (k>=fnPixel3))
                {   
                    LERROR(lmclog,"Imported mode map dimensions don't agree with those in \"" << aFilename <<".\" Double check dim[1,2,3]-max.");
                    return false;
                }   
                std::vector E_input = {Erho,Etheta,Ez};
                fModeMapTE_E[i][j][k] = E_input;
//              printf("read var1 is %g, var2 is %g, E is %g\n", fModeMapTE_E.back()[0], fModeMapTE_E.back()[1], fModeMapTE_E.back()[2]);
            }
        }

        modeMapFile.close();

	//Reset dimensions from import file to actual cavity dimensions in case they don't match up
        MatchCavityDimensions(aParam);

        return true;
    }

    std::vector<double> ModeMapCavity::TE_E(double dim1, double dim2, double dim3, int l, int m, int n, double var1, double var2, double zKass, bool includeOtherPols)
    {
        std::vector<double> TE_E;

        double var3 = zKass + 0.5*dim3;
        std::vector< int > CoordinateIndices = FindClosestCoordinate(var1, var2, var3);
        std::vector< std::vector< int >> TetrahedronVertices = GetVerticesIndices(CoordinateIndices, var1, var2, var3);
        if(CoordinateIndices[0]!=0 or CoordinateIndices[1]!=0 or CoordinateIndices[2]!=0)
        {
            // Found a near neighbor in the uploaded mode map:
            TE_E.push_back( InterpolateField(var1, var2, var3, TetrahedronVertices, 0));  // r
            TE_E.push_back( InterpolateField(var1, var2, var3, TetrahedronVertices, 2));  // z
            TE_E.push_back( InterpolateField(var1, var2, var3, TetrahedronVertices, 1));  // theta

            return TE_E;  // Return the near neighbor.
        }
        else
        {
            std::vector< double > zeroVector(3,0.);
            return zeroVector;
        }
    }

    //These next components need to be expanded once we have real Field Maps to utilize
    std::vector<double> ModeMapCavity::TE_H(double dim1, double dim2, double dim3, int l, int m, int n, double var1, double var2, double zKass, bool includeOtherPols)
    {
//    	LPROG( lmclog, "TE_H is presently not available as a mode map." );
        std::vector<double> TE_H;
        TE_H.push_back(0.);
        return TE_H;
    }

    std::vector<double> ModeMapCavity::TM_E(double dim1, double dim2, double dim3, int l, int m, int n, double var1, double var2, double zKass, bool includeOtherPols)
    {
//      LPROG( lmclog, "TM_E is presently not available as a mode map." );
        std::vector<double> TM_E;
        TM_E.push_back(0.);
        return TM_E;
    }

    std::vector<double> ModeMapCavity::TM_H(double dim1, double dim2, double dim3, int l, int m, int n, double var1, double var2, double zKass, bool includeOtherPols)
    {
//      LPROG( lmclog, "TM_H is presently not available as a mode map." );
        std::vector<double> TM_H;
        TM_H.push_back(0.);
        return TM_H;
    }


    std::vector< int > ModeMapCavity::FindClosestCoordinate(double var1, double var2, double var3)
    {
        //Finds coordinate indices with the floor of the index closest to that input variable for each dimension. Assumes a uniform grid in each of the 3 dimensions.
        std::vector< int > Coordinates(3);
        if( var1<fDim1_min or var1>fDim1_max or var2<fDim2_min or var2>fDim2_max or var3<fDim3_min or var3>fDim3_max )
        {
            Coordinates[0] = 0;
            Coordinates[1] = 0;
            Coordinates[2] = 0;
            return Coordinates;
        }
        else
        {
            Coordinates[0] = (int)((var1 - fDim1_min)/(fDim1_max - fDim1_min)*(fnPixel1-1) + 0.5 - (var1<0)); //  "+ 0.5 - (var1<0)" means casting to int will round to the nearest int rather than truncate towards zero
            Coordinates[1] = (int)((var2 - fDim2_min)/(fDim2_max - fDim2_min)*(fnPixel2-1) + 0.5 - (var2<0));
            Coordinates[2] = (int)((var3 - fDim3_min)/(fDim3_max - fDim3_min)*(fnPixel3-1) + 0.5 - (var3<0));
            return Coordinates;
        }
    }

    std::vector< std::vector< int >> ModeMapCavity::GetVerticesIndices(std::vector<int> ClosestCoordinate, double var1, double var2, double var3)
    {
        //Assumes uniform grid spacing, but gives indicies for the coordinates of the 4 points from fModeMapTE_E that define the smallest tetrahedron enclosing a point in space
        double Closest_v1 = IndexToCoordinate(ClosestCoordinate[0], fDim1_min, fDim1_max, fnPixel1);
        double Closest_v2 = IndexToCoordinate(ClosestCoordinate[1], fDim2_min, fDim2_max, fnPixel2);
        double Closest_v3 = IndexToCoordinate(ClosestCoordinate[2], fDim3_min, fDim3_max, fnPixel3);
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

    double ModeMapCavity::InterpolateField(double var1, double var2, double var3, std::vector< std::vector<int>> TetrahedronVertices, int component)
    {
        //Does linear interpolation of a field at a point within a tetrahedron with known field at the vertices
        double x0, x1, x2, x3, y0, y1, y2, y3, z0, z1, z2, z3;
        x0 = IndexToCoordinate(TetrahedronVertices[0][0], fDim1_min, fDim1_max, fnPixel1);
        x1 = IndexToCoordinate(TetrahedronVertices[1][0], fDim1_min, fDim1_max, fnPixel1);
        x2 = IndexToCoordinate(TetrahedronVertices[2][0], fDim1_min, fDim1_max, fnPixel1);
        x3 = IndexToCoordinate(TetrahedronVertices[3][0], fDim1_min, fDim1_max, fnPixel1);
        y0 = IndexToCoordinate(TetrahedronVertices[0][1], fDim2_min, fDim2_max, fnPixel2);
        y1 = IndexToCoordinate(TetrahedronVertices[1][1], fDim2_min, fDim2_max, fnPixel2);
        y2 = IndexToCoordinate(TetrahedronVertices[2][1], fDim2_min, fDim2_max, fnPixel2);
        y3 = IndexToCoordinate(TetrahedronVertices[3][1], fDim2_min, fDim2_max, fnPixel2);
        z0 = IndexToCoordinate(TetrahedronVertices[0][2], fDim3_min, fDim3_max, fnPixel3);
        z1 = IndexToCoordinate(TetrahedronVertices[1][2], fDim3_min, fDim3_max, fnPixel3);
        z2 = IndexToCoordinate(TetrahedronVertices[2][2], fDim3_min, fDim3_max, fnPixel3);
        z3 = IndexToCoordinate(TetrahedronVertices[3][2], fDim3_min, fDim3_max, fnPixel3);

        for(int i=0; i<4; i++) //check if any indices are outside the grid of points used for interpolation
        {
            if(TetrahedronVertices[i][0]>=fnPixel1 or TetrahedronVertices[i][2]>=fnPixel3) return 0.; //if r or z are outside of the grid, set field to zero
            if(TetrahedronVertices[i][0]<0 or TetrahedronVertices[i][2]<0) return 0.;
            if(TetrahedronVertices[i][1]>=fnPixel2) TetrahedronVertices[i][1] -= fnPixel2; //if theta is outside the the grid size wrap around back to 0 for periodicity
            if(TetrahedronVertices[i][1]<0) TetrahedronVertices[i][1] += fnPixel2;
        }

        Eigen::MatrixXd m {
            {1., x0, y0, z0},
            {1., x1, y1, z1},
            {1., x2, y2, z2},
            {1., x3, y3, z3}};

        Eigen::VectorXd v {{ fModeMapTE_E[TetrahedronVertices[0][0]][TetrahedronVertices[0][1]][TetrahedronVertices[0][2]][component], fModeMapTE_E[TetrahedronVertices[1][0]][TetrahedronVertices[1][1]][TetrahedronVertices[2][2]][component], fModeMapTE_E[TetrahedronVertices[2][0]][TetrahedronVertices[2][1]][TetrahedronVertices[2][2]][component], fModeMapTE_E[TetrahedronVertices[3][0]][TetrahedronVertices[3][1]][TetrahedronVertices[3][2]][component]}};

        Eigen::VectorXd Coef = m.inverse() * v;

        Eigen::VectorXd pos {{ 1., var1, var2, var3}};

        double interpolated_value = pos.dot(Coef);
        return interpolated_value;
    }

    double ModeMapCavity::IndexToCoordinate(int index, double min, double max, int nPixels)
    {
        return (double)index*(max - min)/((double)(nPixels-1)); //nPixels +1 since range in text file goes from bin midpoint to midpoint
    }

    void ModeMapCavity::MatchCavityDimensions(const scarab::param_node& aParam)
    {

        if( aParam.has( "cavity-radius" ) ) 
        {   
            fDim1_max = aParam["cavity-radius"]().as_double() ;
        }   
        else if( aParam.has( "cavity-x" ) )
        {
      	    fDim1_max = aParam["cavity-x"]().as_double();
        }

        if( aParam.has( "cavity-y" ) )
        {   
            fDim2_max = aParam["cavity-y"]().as_double() ;
        }  	
        if( aParam.has( "cavity-length" ) ) 
        {   
            fDim3_max = aParam["cavity-length"]().as_double() ;
        }  
    }

} /* namespace locust */

