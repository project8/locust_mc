/*
 * LMCFileWriter.cc
 *
 *  Created on: Mar 21, 2020
 *      Author: pslocum
 */

#include "LMCFileWriter.hh"

namespace locust
{
	FileWriter::FileWriter() {}
    FileWriter::~FileWriter() {}

    bool FileWriter::Configure( const scarab::param_node& aParam )
     {
          return true;
     }

    /*
    double FileWriter::GetTestVar()
    {
    	return fTestVar;
    }

    void FileWriter::SetTestVar(double aValue)
    {
    	fTestVar = aValue;
    }
*/


} /* namespace locust */

