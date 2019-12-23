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

    std::vector<std::deque<double>> empty;
    empty.resize(buffer.size());

    for (unsigned index = 0; index < buffer.size(); index ++ )
    {
    	for (std::deque<double>::iterator itdeque = buffer[index].begin(); itdeque != buffer[index].end(); ++itdeque)
    	{
    		empty[index].push_back(buffer[index].front());
    		buffer[index].pop_front();
    	}
    }

    std::swap( buffer, empty );
    return buffer;
    }


  std::vector<std::deque<unsigned>> FieldBuffer::CleanupBuffer(std::vector<std::deque<unsigned>> buffer)
    {

	std::vector<std::deque<unsigned>> empty;
	empty.resize(buffer.size());

	for (unsigned index = 0; index < buffer.size(); index ++ )
	{
	    for (std::deque<unsigned>::iterator itdeque = buffer[index].begin(); itdeque != buffer[index].end(); ++itdeque)
	    {
	    	empty[index].push_back(buffer[index].front());
	    	buffer[index].pop_front();
	    }
	}

	std::swap( buffer, empty );
    return buffer;

    }



} /* namespace locust */
