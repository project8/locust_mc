/*
 * LMCFileWriter.cc
 *
 *  Created on: Mar 21, 2020
 *      Author: pslocum
 */

#include "LMCFileWriter.hh"

namespace locust
{
    FileWriter::FileWriter()
	{
	}
 
    FileWriter::~FileWriter() {}

    bool FileWriter::Configure( const scarab::param_node& aParam )
    {
        return true;
    }

    void FileWriter::SetFilename(std::string aFilename)
    {
    	fRoot_filename = aFilename;
    }

    std::string FileWriter::GetFilename()
    {
    	return fRoot_filename;
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

