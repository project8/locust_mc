#include <map>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
using namespace std;
#include <stdio.h>
#include <string.h>
#include <yajl/yajl_gen.h>
#include <yajl/yajl_tree.h>


class TransferFunction
{
public:
    TransferFunction();
    TransferFunction(const TransferFunction &f);
    TransferFunction(double fstart,double fspan,unsigned int sz);
    ~TransferFunction();
    TransferFunction &operator=(const TransferFunction &f);
    TransferFunction interpolateTo(double fstart,double fspan,unsigned int sz);
    float interpolatePoint(double freq);
    void resize(unsigned int sz); //note, destroys data

    void loadFromFile(string fname);
    void saveToJson(yajl_gen node);
    void loadFromJson(yajl_val node);

    float *data;
    unsigned int size;
    double frequency_start;
    double frequency_span;

};

class ReceiverTransferFunctions
{
public:
    TransferFunction high_frequency_stage;
    map<float,TransferFunction> low_frequency_stage; //mixing freq, function
    double hf_lo_freq; //in MHz

    TransferFunction getTransferFunction(double lo_freq,double fspan,unsigned int sz);

    void load_from_json(string fname);
    void load_lf_from_filelist(string fname);
    void save_to_file(string fname);
};
