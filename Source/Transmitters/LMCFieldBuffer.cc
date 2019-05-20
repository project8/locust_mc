/*
 * LMCFieldBuffer.cc
 *
 *  Created on: May 19, 2019
 *      Author: pslocum
 */

#include "LMCFieldBuffer.hh"


namespace locust
{

    FieldBuffer::FieldBuffer():
    buffer(1)
    {
    }

    FieldBuffer::~FieldBuffer()
    {
    }


    std::vector<std::deque<double>> FieldBuffer::InitializeBuffer(int nchannels, int buffersize)
    {
    std::vector<std::deque<double>> buffer;
    buffer.resize(nchannels);
    for (int i=0; i<nchannels; i++)
    {
    for (int j=0; j<buffersize; j++)
      {
      buffer[i].emplace(buffer[i].begin()+j, 0.);
      }
    return buffer;
    }
    }




} /* namespace locust */
