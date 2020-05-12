/*
 * LMCFieldParameterGenerator.hh
 *
 *  Created on: March 04, 2020
 *      Author: P. T. Surukuchi
 */

#ifndef LMCFIELDPARAMETERGENERATOR_HH_
#define LMCFIELDPARAMETERGENERATOR_HH_

#include "LMCTransmitterInterfaceGenerator.hh"

#include "LMCEventHold.hh"
#include "LMCRunKassiopeia.hh"

#include "logger.hh"
#include <thread>
#include <algorithm>

#include <iostream>
#include <fstream>

namespace locust
{

    /*!
      @class FieldParameterGenerator
      @author P. T. Surukuchi

      @brief Generate signal in free space(without wave guide) for phase III and estimate field parameters (Poynting vector etc) 

      @details
      Operates in time space
      Currently two possible ways to define the field points
      1. Use a text file 
      2. Predefined 

      When using predefined points, the options currently are:
      1. A sphere(1) centered at the origin
      2. A cylinder(2) centered at the origin

      Configuration name: "field-parameter"
      - "use-text-file": bool -- Flag to check if a text file will be input to define the field points
      - "predefined-geometry": string -- Integer value for defining the predefined geometry to be used to generate field points
      - "file-name": string -- input textfile with the points in space in cartesian coordinate system where the field points need to be defined 

      Available configuration options:

*/

    class FieldParameterGenerator : public TransmitterInterfaceGenerator 
    {
        public:

            FieldParameterGenerator( const std::string& aName = "field-parameter" );
            virtual ~FieldParameterGenerator();

            bool Configure( const scarab::param_node& aNode );

        private:
            bool fUseTextFile;
            int fPredefinedGeometry;
            double fRadius;
            double fLength;
            double fMinFieldPoints;
            std::string fTextFileName;
            std::vector<LMCThreeVector> fFieldPoints;
            std::vector<std::vector<double>> fTotalField;

            bool ends_with(const std::string &str, const std::string &suffix);
            bool DoGenerate(Signal* aSignal);
            void InitializeFieldPoints();
            double GenerateFieldPoints();
            void DriveAntenna(FILE *fp, int PreEventCounter, unsigned index, Signal* aSignal, int nfilterbins, double dtfilter);
    };

} /* namespace locust */

#endif /* LMCFIELDPARAMETERGENERATOR_HH_ */
