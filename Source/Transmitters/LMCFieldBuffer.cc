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


    std::vector<std::deque<double>> FieldBuffer::InitializeBuffer(int nchannels, int npatches, int buffersize)
    {
    std::vector<std::deque<double>> buffer;
    buffer.resize(nchannels*npatches);
    for (int i=0; i<nchannels; i++)
      {
      for (int j=0; j<npatches; j++)
        {
        for (int k=0; k<buffersize; k++)
          {
          buffer[i*npatches+j].emplace(buffer[i*npatches+j].begin()+k, 0.);
          }
        }
      }
    return buffer;
    }




} /* namespace locust */
