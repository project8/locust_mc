
#include "Monarch.hpp"
#include <stdlib.h>
#include <string.h>
#include <iostream>
using namespace std;


int main(int argc,char *argv[])
{
    char *eggname=argv[1];
    double histogram[256];
    for(int i=0;i<256;i++) {
	histogram[i]=0;
    }
    const Monarch *egg=Monarch::OpenForReading(std::string(eggname));
    egg->ReadHeader();
    const MonarchHeader *eggheader=egg->GetHeader();
    const MonarchRecord *event;
//    sampling_rate_mhz=eggheader->GetAcqRate();
    while(egg->ReadRecord()) {
	event=egg->GetRecordOne();
      	for(int i=0;i<eggheader->GetRecordSize();i++) {
	    histogram[event->fDataPtr[i]]+=1;
	}
    }
    for(int i=0;i<256;i++) {
	cout << i << " " << histogram[i] << endl;
    }
}
