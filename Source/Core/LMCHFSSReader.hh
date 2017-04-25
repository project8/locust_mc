#ifndef LMCHFSSREADER_HH_
#define LMCHFSSREADER_HH_

#include "LMCException.hh"

#include <string>
#include <vector>
#include <array>
#include <sstream>

namespace locust
{

    class HFSSReader 
    {
        public:
            HFSSReader();
            void ParseANDFile(std::string fAND_filename);
            std::vector<std::array<double, 3> > GetSurfacePoints();
            std::vector<std::array<double, 3> > GeneratePlane(std::array<double, 2> GeometryScale, int nResolution);
            std::vector<std::array<double, 3> > RotateShift(std::vector<std::array<double, 3> > rPointVector, std::array<double, 3> tNormal, std::array<double, 3> rCenter);
            std::vector<double> GetFrequencies();
            std::string GetNFDFilename();

        private:
            //std::string fAND_filename;
            std::vector<std::array<double, 3> > rSurfacePoints; //vector of 3D positions of points at which the fields are calculated at. Points form one of surrfaces below
            std::array<double, 3> GeometryCenter; 
            std::array<double, 3> GeometryScale; 
            std::array<double, 3> GeometryAxis; 

            std::string fNFD_filename;
            std::vector<double> NFDFrequencies;

            double ParseUnits(std::string TextInput); //Gets numerical value for units
            void StringClean(std::string &InputString); //Removes commas/ unnecessary text from line
            void ArrayParse(std::string InputString, std::array<double, 3> (&X) ); //Parses lists of numbers: "a, b, c"

            std::vector<std::array<double, 3> > GenerateSphere(double Radius, int nResolution);
            std::vector<std::array<double, 3> > GenerateBox(std::array<double, 3> GeometryScale, int nResolution);
            std::vector<std::array<double, 3> > GenerateCylinder(std::array<double, 2> GeometryScale, int nResolution);


    };

}

#endif
