/*
 * LMCHFSSReader.cc
 *
 *  Created on: Apr 16, 2017
 *      Author: nbuzinsky
 */

#ifndef PI
#define PI 3.1415926
#endif

#include "LMCHFSSReader.hh"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <regex>
#include <math.h>

#include "logger.hh"

namespace locust
{
    LOGGER( lmclog, "HFSSReader" );

    HFSSReader::HFSSReader():
        fNFD_filename("freefield.nfd"),
        GeometryCenter({0.,0.,0.}),
        GeometryScale({1.,1.,1.}),
        GeometryAxis({0.,0.,1.})
    {
    }

    std::vector<double> HFSSReader::GetFrequencies(){return NFDFrequencies;}

    std::string HFSSReader::GetNFDFilename(){return fNFD_filename;}

    void HFSSReader::ParseANDFile(std::string fAND_filename)
    {
        std::ifstream fANDInput;
        fANDInput.open(fAND_filename,std::ios::in);
        std::string InputLine;

        std::regex ParseExpressionHeader;
        std::regex ParseExpressionLine;
        std::smatch ParseMatch;

        ParseExpressionHeader="\\$begin .NearFieldHeader";

        //Loop over empty lines
        do{
            getline(fANDInput,InputLine);
        }while(!std::regex_search(InputLine,ParseMatch,ParseExpressionHeader));

        //Get Data from NearFieldHeader
        ParseExpressionHeader="\\$end .NearFieldHeader";
        ParseExpressionLine="[:alpha:]+.*=.*";

        int GeometryIndex=0; //0-Sphere 1-Box 2-Cylinder 3-Plane

        do{
            getline(fANDInput,InputLine);
            //std::cout<<InputLine<<std::endl;

            if(InputLine.find("type")!=std::string::npos)
            {
                if(InputLine.find("nfd")==std::string::npos)
                {
                    ERROR( lmclog, "Must produce NFD file if AND input is specified!" );
                }
            }
            else if(InputLine.find("fields")!=std::string::npos)
            {
                if(InputLine.find("EH")==std::string::npos)
                {
                    ERROR( lmclog, "Program only supports EH Field Mode!" );
                }
            }
            else if(InputLine.find("geometry")!=std::string::npos)
            {
                if(InputLine.find("sphere")!=std::string::npos) GeometryIndex=0;
                else if(InputLine.find("box")!=std::string::npos) GeometryIndex=1;
                else if(InputLine.find("cylinder")!=std::string::npos) GeometryIndex=2;
                else if(InputLine.find("plane")!=std::string::npos) GeometryIndex=3;
                else
                {
                    ERROR( lmclog, "Program only supports sphere, box, cylinder or plane!" );
                }
            }
            else if(InputLine.find("center")!=std::string::npos)
            {
                ArrayParse(InputLine,GeometryCenter);
            }
            else if(InputLine.find("size")!=std::string::npos)
            {
                ArrayParse(InputLine,GeometryScale);
            }

            else if(InputLine.find("axis")!=std::string::npos)
            {
                ArrayParse(InputLine,GeometryAxis);
            }

            else if(InputLine.find("normal")!=std::string::npos)
            {
                ArrayParse(InputLine,GeometryAxis);
            }

            else if(InputLine.find("radius")!=std::string::npos)
            {
                StringClean(InputLine);
                double fRadius=ParseUnits(InputLine);
                GeometryScale[0]=fRadius;
            }

            else if(InputLine.find("height")!=std::string::npos)
            {
                StringClean(InputLine);
                double fHeight=ParseUnits(InputLine);
                GeometryScale[1]=fHeight;
            }

            else if(InputLine.find("fsweep")!=std::string::npos)
            {
                StringClean(InputLine);
                std::stringstream InputStream(InputLine);
                std::string InputBuffer;
                while(InputStream >> InputBuffer)
                    NFDFrequencies.push_back(ParseUnits(InputBuffer));
            }

            else if(InputLine.find("fcoords")!=std::string::npos)
            {
                if(InputLine.find("cartesian")==std::string::npos)
                    ERROR( lmclog, "Program only supports cartesian coordinates for now!" );
            }


        }while(!std::regex_search(InputLine,ParseMatch,ParseExpressionHeader));


        //Loop over empty lines
        ParseExpressionHeader="\\$begin .NearFieldData";
        do{
            getline(fANDInput,InputLine);
            //std::cout<<InputLine<<std::endl;
        }while(!std::regex_search(InputLine,ParseMatch,ParseExpressionHeader));

        
        //Get NFD File Name
        ParseExpressionHeader="\\$end .NearFieldData";
        do{
            getline(fANDInput,InputLine);
            if(InputLine.find(".nfd")!=std::string::npos)break;
        }while(!std::regex_search(InputLine,ParseMatch,ParseExpressionHeader));

        ParseExpressionLine=".*nfd";
        std::regex_search(InputLine,ParseMatch,ParseExpressionLine);
        std::string fNFD_filename=ParseMatch.str();
        //Get out Filename
        int tIndex=0;
         while((tIndex=fNFD_filename.find("\""))!=std::string::npos)
         {
             fNFD_filename.erase(0,tIndex+1);
         }
         //cout<<fNFD_filename<<std::endl;

        fANDInput.close();

        //Generate Geometries: Increase Final Argument for finer discretization 
        if(GeometryIndex==0)rSurfacePoints=GenerateSphere(GeometryScale[0],10);
        else if(GeometryIndex==1)rSurfacePoints=GenerateBox(GeometryScale,13);
        else if(GeometryIndex==2)rSurfacePoints=GenerateCylinder({GeometryScale[0],GeometryScale[1]},10);
        else if(GeometryIndex==3)rSurfacePoints=GeneratePlane({GeometryScale[0],GeometryScale[1]},10);
     
        //Print DEBUG information. Make sure that everything looks correct
        DEBUG( lmclog, "Configuring HFSS Compatible Output");

        if(GeometryIndex==0)
        {
            DEBUG( lmclog, "Sphere generated with "<<rSurfacePoints.size()<<" surface points");
        }
        else if(GeometryIndex==1)
        {
            DEBUG( lmclog, "Box generated with "<<rSurfacePoints.size()<<" surface points");
        }
        else if(GeometryIndex==2)
        {
            DEBUG( lmclog, "Cylinder generated with "<<rSurfacePoints.size()<<" surface points");
        }
        else if(GeometryIndex==3)
        {
            DEBUG( lmclog, "Plane generated with "<<rSurfacePoints.size()<<" surface points");
        }

        DEBUG( lmclog, "Geometry Center: "<<GeometryCenter[0]<<" "<<GeometryCenter[1]<<" "<<GeometryCenter[2] );
        DEBUG( lmclog, "Geometry Scale: "<<GeometryScale[0]<<" "<<GeometryScale[1]<<" "<<GeometryScale[2] );

        DEBUG( lmclog, "Chosen HFSS Frequencies:");
        for(int i=0;i<NFDFrequencies.size();i++) DEBUG( lmclog, NFDFrequencies[i]);

        return;
    }

    std::vector<std::array<double, 3> > HFSSReader::GetSurfacePoints()
    {
        return rSurfacePoints;
    }

    double HFSSReader::ParseUnits(std::string TextInput)
    {
        double RealValue=0.;
        double Scale=1.;
        int tIndex;

        if((tIndex=TextInput.find("mm"))!=std::string::npos)Scale=1e-3;
        else if((tIndex=TextInput.find("nm"))!=std::string::npos)Scale=1e-9;
        else if((tIndex=TextInput.find("meter"))!=std::string::npos)Scale=1.;
        else if((tIndex=TextInput.find("cm"))!=std::string::npos)Scale=1e-2;
        else if((tIndex=TextInput.find("ft"))!=std::string::npos)Scale=0.3048;
        else if((tIndex=TextInput.find("m"))!=std::string::npos)Scale=1.;
        else if((tIndex=TextInput.find("in"))!=std::string::npos)Scale=0.0254;
        else if((tIndex=TextInput.find("kHz"))!=std::string::npos)Scale=1e3;
        else if((tIndex=TextInput.find("MHz"))!=std::string::npos)Scale=1e6;
        else if((tIndex=TextInput.find("GHz"))!=std::string::npos)Scale=1e9;
        else if((tIndex=TextInput.find("THz"))!=std::string::npos)Scale=1e12;
        else if((tIndex=TextInput.find("Hz"))!=std::string::npos)Scale=1.;

        if(tIndex!=std::string::npos)
        {
            TextInput.resize(tIndex);
        }

        RealValue=std::stod(TextInput)*Scale;
        return RealValue;
    }

    void HFSSReader::StringClean(std::string &InputString)
    {
        int StringIndex=InputString.find("'"); 
        InputString.erase(0,StringIndex+1);
        InputString.erase(InputString.size()-1);
        InputString.erase(std::remove(InputString.begin(), InputString.end(), ','), InputString.end());
    }

    void HFSSReader::ArrayParse(std::string InputString, std::array<double, 3> (&X) )
    {
        StringClean(InputString);
        std::stringstream InputStream(InputString);
        std::string InputBuffer;
        for(int i=0;i<3;i++)
        {
            InputStream >> InputBuffer;
            if(!InputBuffer.size())break;
            X[i]=ParseUnits(InputBuffer);
        }
    }


    std::vector<std::array<double, 3> > HFSSReader::GeneratePlane(std::array<double, 2> GeometryScale, int nResolution)
    {
        std::vector<std::array<double, 3> > rPointVector;

        double SizeRatio[2]={1.,1.};
        int MinIndex=(GeometryScale[0] < GeometryScale[1]) ? 0 : 1;
        SizeRatio[(MinIndex+1)%2]=GeometryScale[(MinIndex+1)%2]/GeometryScale[MinIndex];

        int N[2]; double dx[2];
        for(int i=0;i<2;i++)
        {
            N[i]=int(double(nResolution)*SizeRatio[i])+1;
            dx[i]=double(GeometryScale[i])/double(N[i]-1);
        }

        std::array<double, 3> rPointVectorBuffer;
        rPointVectorBuffer[2]=0.;
        for(int i=0;i<N[0];i++)
        {
            rPointVectorBuffer[0]=dx[0]*double(i);
            for(int j=0;j<N[1];j++)
            {
                rPointVectorBuffer[1]=dx[1]*double(j);
                rPointVector.push_back(rPointVectorBuffer);
            }
        }
        rPointVector=RotateShift(rPointVector,{0.,0.,1.},{-dx[0]*double(N[0]-1.)/2.,-dx[1]*double(N[1]-1.)/2.,0.});

        return rPointVector;
    }

    std::vector<std::array<double, 3> > HFSSReader::GenerateBox(std::array<double, 3> GeometryScale, int nResolution)
    {
        std::vector<std::array<double, 3> > rPointVector;

        double SizeRatio[3]={1.,1.,1.};
        int MinIndex=0.;
        for(int i=1;i<3;i++)
            if(GeometryScale[i]<GeometryScale[MinIndex]) MinIndex=i;

        SizeRatio[(MinIndex+1)%3]=GeometryScale[(MinIndex+1)%3]/GeometryScale[MinIndex];
        SizeRatio[(MinIndex+2)%3]=GeometryScale[(MinIndex+2)%3]/GeometryScale[MinIndex];

        int N[3]; double dx[3];
        for(int i=0;i<3;i++)
        {
            N[i]=int(double(nResolution)*SizeRatio[i])+1;
            dx[i]=double(GeometryScale[i])/double(N[i]-1);
        }

        std::array<double, 3> rPointVectorBuffer;
        for(int i=0;i<N[0];i++)
        {
            rPointVectorBuffer[0]=dx[0]*(double(i)-double(N[0]-1.)/2.);
            for(int j=0;j<N[1];j++)
            {
                rPointVectorBuffer[1]=dx[1]*(double(j)-double(N[1]-1.)/2.);
                for(int k=0;k<N[2];k++)
                {
                    rPointVectorBuffer[2]=dx[2]*(double(k)-double(N[2]-1.)/2.);
                    if(abs(rPointVectorBuffer[0])==GeometryScale[0]/2. || abs(rPointVectorBuffer[1])==GeometryScale[1]/2. || abs(rPointVectorBuffer[2])==GeometryScale[2]/2.)
                        rPointVector.push_back(rPointVectorBuffer);
                }
            }
        }

        return rPointVector;

    }

    std::vector<std::array<double, 3> > HFSSReader::GenerateSphere(double Radius, int nResolution)
    {
        std::vector<std::array<double, 3> > rPointVector;
        int nRings=2.*nResolution+3;
        double dz=2./double(nRings-1.);
        double zRings[nRings];
        double rRings[nRings];
        double RingCount[nRings];

        for(int i=0;i<nRings;i++)
        {
            zRings[i]=-1.+dz*double(i);
            rRings[i]=sqrt(1.-pow(zRings[i],2.));
            RingCount[i]=round(2.*PI*rRings[i])/dz;
        }
        double Phi=0.;

        std::array<double, 3> rPointVectorBuffer;
        for(int i=0;i<nRings;i++)
        {
            Phi=0.;
            rPointVectorBuffer[2]=Radius*zRings[i];
            double dPhi=2.*PI/double(RingCount[i]);
            //if(i%2)Phi+=dPhi/3.;
            do{
                rPointVectorBuffer[0]=Radius*rRings[i]*cos(Phi);
                rPointVectorBuffer[1]=Radius*rRings[i]*sin(Phi);
                rPointVector.push_back(rPointVectorBuffer);
                Phi+=dPhi;
            }while(Phi<2.*PI);
        }

        return rPointVector;
    }

    std::vector<std::array<double, 3> > HFSSReader::GenerateCylinder(std::array<double, 2> GeometryScale, int nResolution)
    {
        std::vector<std::array<double, 3> > rPointVector;

        double SizeRatio[2]={1.,1.};
        int MinIndex=(GeometryScale[0] < GeometryScale[1]) ? 0 : 1;
        SizeRatio[(MinIndex+1)%2]=GeometryScale[(MinIndex+1)%2]/GeometryScale[MinIndex];

        int N[2]; double dz[2];
        for(int i=0;i<2;i++)
        {
            N[i]=int(double(nResolution)*SizeRatio[i])+1;
            dz[i]=double(GeometryScale[i])/double(N[i]-1);

        }

        double Radius=GeometryScale[0];
        int nRings=N[1];
        int RingCount=int(2.*PI*Radius/dz[1]);
        double dPhi=2.*PI/double(RingCount);
        double Phi=0.; double Z=0.;
        //Fill in outer surface of cylinder
        std::array<double, 3> rPointVectorBuffer;
        for(int i=0;i<nRings;i++)
        {
            Phi=0.;
            rPointVectorBuffer[2]=Z;
            for(int j=0;j<RingCount;j++)
            {
                rPointVectorBuffer[0]=Radius*cos(Phi);
                rPointVectorBuffer[1]=Radius*sin(Phi);
                rPointVector.push_back(rPointVectorBuffer);
                Phi+=dPhi;
            }
            Z+=dz[1];
        }
        //Fill in end caps
        std::vector<std::array<double, 3> > rPointVectorEnd;
        rPointVectorBuffer[0]=0.; rPointVectorBuffer[1]=0.; rPointVectorBuffer[2]=0.;
        rPointVectorEnd.push_back(rPointVectorBuffer);
        for(int i=1;i<N[0]-1;i++)
        {
            RingCount=round(2.*PI*double(i));
            Phi=0.;
            dPhi=2.*PI/double(RingCount);
            while(Phi<2.*PI)
            {
                rPointVectorBuffer[0]=dz[0]*double(i)*cos(Phi);
                rPointVectorBuffer[1]=dz[0]*double(i)*sin(Phi);
                rPointVectorEnd.push_back(rPointVectorBuffer);
                Phi+=dPhi;

            }
        }
        rPointVector.insert(rPointVector.end(), rPointVectorEnd.begin(), rPointVectorEnd.end());
        rPointVectorEnd=RotateShift(rPointVectorEnd,{0.,0.,1.},{0.,0.,GeometryScale[1]});
        rPointVector.insert(rPointVector.end(), rPointVectorEnd.begin(), rPointVectorEnd.end());
        rPointVector=RotateShift(rPointVector,{0.,0.,1.},{0.,0.,-GeometryScale[1]/2.});

        return rPointVector;

    }

    std::vector<std::array<double, 3> > HFSSReader::RotateShift(std::vector<std::array<double, 3> > rPointVector, std::array<double, 3> tNormal, std::array<double, 3> rCenter)
    {
        double tNormalization=0.;
        for(int i=0;i<3;i++)tNormalization+=pow(tNormal[i],2.);
        tNormalization=sqrt(tNormalization);
        if(tNormalization==0.)
        {
            tNormalization=1.;
            tNormal[2]=1.;
        }
        for(int i=0;i<3;i++)tNormal[i]/=tNormalization;


        //Deal with infinities
        double Theta=acos(tNormal[2]);

        //If Both XY normals are 0, do not do phi rotation
        double Phi=0.;
        if(tNormal[0])Phi=atan(tNormal[1]/tNormal[0]);
        else if(tNormal[1])Phi=PI/2.;


        ///Perform Rotation on surface
        double SurfaceRotation[3][3]={{cos(Theta)*cos(Phi),-sin(Phi),cos(Phi)*sin(Theta)},{cos(Theta)*sin(Phi),cos(Phi),sin(Phi)*sin(Theta)},{-sin(Theta),0.,cos(Theta)}};

        for(unsigned ix=0;ix<rPointVector.size();ix++)
        {
            double tmpPos[3]={};
            for(int i=0;i<3;i++)
            {
                for(int j=0;j<3;j++)
                {
                    tmpPos[i]+=SurfaceRotation[i][j]*rPointVector[ix][j];
                }
            }
            //Shift Surface to desired location
            for(int i=0;i<3;i++)
            {
                rPointVector[ix][i]=tmpPos[i]+rCenter[i];
            }
        }

        return rPointVector;
    }

} /* namespace locust */
