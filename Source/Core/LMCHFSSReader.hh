#ifndef LMCHFSSREADER_HH_
#define LMCHFSSREADER_HH_

#include "LMCThreeVector.hh"
#include "LMCException.hh"

#include <string>
#include <vector>
#include <array>
#include <sstream>

namespace locust
{
    /*!
         @class HFSSReader
         @author N. Buzinsky
         @brief Parses the HFSS .and file
         @details Peak phasor quantities (specific frequencies) at 3D positions of points at which the fields are calculated. Points on surface of sphere, box, or cylinder
    */
    
    class HFSSReader 
    {
        public:
            HFSSReader();
            void ParseANDFile(std::string fAND_filename);
            std::vector< LMCThreeVector> GetSurfacePoints();
            std::vector< LMCThreeVector> GeneratePlane(std::array<double, 2> GeometryScale, int nResolution);
            std::vector< LMCThreeVector> RotateShift(std::vector<LMCThreeVector> rPointVector, LMCThreeVector tNormal, LMCThreeVector rCenter);
            std::vector<double> GetFrequencies();
            std::string GetNFDFilename();

        private:
            //std::string fAND_filename;
            std::vector< LMCThreeVector> rSurfacePoints; //vector of 3D positions of points at which the fields are calculated at. Points form one of surfaces below
            LMCThreeVector GeometryCenter; 
            LMCThreeVector GeometryScale; 
            LMCThreeVector GeometryAxis; 

            std::string fNFD_filename;
            std::vector<double> NFDFrequencies;

            double ParseUnits(std::string TextInput); //Gets numerical value for units
            void StringClean(std::string &InputString); //Removes commas/ unnecessary text from line
            void ArrayParse(std::string InputString, LMCThreeVector (&X) ); //Parses lists of numbers: "a, b, c"

            std::vector< LMCThreeVector> GenerateSphere(double Radius, int nResolution);
            std::vector< LMCThreeVector> GenerateBox(LMCThreeVector GeometryScale, int nResolution);
            std::vector< LMCThreeVector> GenerateCylinder(std::array<double, 2> GeometryScale, int nResolution);


    };

}

#endif
