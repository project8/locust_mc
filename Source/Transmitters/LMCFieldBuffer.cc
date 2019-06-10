/*
 * LMCFieldBuffer.cc
 *
 *  Created on: May 19, 2019
 *      Author: pslocum
 */

#include "LMCFieldBuffer.hh"


namespace locust
{

    FieldBuffer::FieldBuffer()
    {
    }

    FieldBuffer::~FieldBuffer()
    {
    }


    std::vector<std::deque<double>> FieldBuffer::InitializeBuffer(int nchannels, int npatches, int buffersize)
    {
    std::vector<std::deque<double>> buffer;
    buffer.resize(nchannels*npatches);
    for (int i=0; i<nchannels; i++)
      {
      for (int j=0; j<npatches; j++)
        {
    	for (int k=0; k<buffersize; ++k)
    		buffer[i*npatches+j].push_back(0.0);
        }
      }
    return buffer;
    }


    std::vector<std::deque<unsigned>> FieldBuffer::InitializeUnsignedBuffer(int nchannels, int npatches, int buffersize)
    {
    std::vector<std::deque<unsigned>> buffer;
    buffer.resize(nchannels*npatches);
    for (int i=0; i<nchannels; i++)
      {
      for (int j=0; j<npatches; j++)
        {
    	for (int k=0; k<buffersize; ++k)
    	{
    		buffer[i*npatches+j].push_back(0);
    	}
        }
      }
    return buffer;
    }





  std::vector<std::deque<double>> FieldBuffer::CleanupBuffer(std::vector<std::deque<double>> buffer)
    {

	for (unsigned i=0; i<buffer.size(); i++)
	{
		buffer[i].clear();
		buffer[i].shrink_to_fit();
	}
    buffer.clear();
    buffer.shrink_to_fit();
    return buffer;
    }


  std::vector<std::deque<unsigned>> FieldBuffer::CleanupBuffer(std::vector<std::deque<unsigned>> buffer)
    {
    for (unsigned i=0; i<buffer.size(); i++)
	  {
	  buffer[i].clear();
      buffer[i].shrink_to_fit();
      }
    buffer.clear();
    buffer.shrink_to_fit();
    return buffer;
    }



} /* namespace locust */
