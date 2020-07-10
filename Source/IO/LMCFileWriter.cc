/*
 * LMCFileWriter.cc
 *
 *  Created on: Mar 21, 2020
 *      Author: pslocum
 */

#include "LMCFileWriter.hh"

namespace locust
{
    FileWriter::FileWriter():
	fFile ( 0 ),
    	fRoot_filename( "LocustEvent.root" )
	{
	}
 
    FileWriter::~FileWriter() {}

    bool FileWriter::Configure( const scarab::param_node& aParam )
    {

    	if( aParam.has( "roothisto-filename" ) )
        {
            fRoot_filename = aParam["roothisto-filename"]().as_string();
        }

        return true;
    }

    void FileWriter::OpenFile(std::string aFileFlag)
    {
        fFile = new TFile(fRoot_filename.c_str(), aFileFlag.c_str());
    }

    void FileWriter::CloseFile()
    {
    	fFile->Close();
    }



} /* namespace locust */

