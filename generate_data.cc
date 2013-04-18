#include <complex.h>
#include "Monarch.hpp"
#include <math.h>
#include <fftw3.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <yajl/yajl_gen.h>
using namespace std;

//program constants
//unsigned char mcname[]="scarab_mc";
char mcname[]="scarab_mc";
//ipunsigned char token_mc_id[]="mc_id";

//physical constants
float kB=1.3806e-23;  // Boltzman constant W/K
float c=3e10; //speed of light in cm/s
double ecyclo=1.758820e11; //electron cyclotron frequency in rad/s
double emass=510998.9; //electron mass in eV

// --<A|---<B>-----<C>----|D>-----
// ------------lwg--------------dl

float c_wg=c*sqrt(1-pow(c/(2.0*0.42*2.54*26e9),2.0)); //speed of light in waveguide in cm/s
float c_cable=0.6*c; //speed of light in cables in cm/s (this is a guess)
float lwg=20; //waveguide length in cm (TODO measure this)
float dl=5; //extra length on channel 2 cable in cm
int record_size=4194304; //number of samples in record time series
int fft_size=record_size/2+1;
float digitizer_fullscale=0.5; //in volts
//float amp_temp=0;
float amp_temp=30;
//float uncor_temp=0;
float uncor_temp=30;
float T_A=uncor_temp; //temperature of noise source A
float T_B=amp_temp; //temperature of noise source B
float T_C=amp_temp; //temperature of noise source C
float T_D=uncor_temp; //temperature of noise source D
float BField=0.90; //in tesla
fftwf_complex *Anoise_f; //A noise in frequency space
fftwf_complex *Bnoise_f; //B noise in frequency space
fftwf_complex *Cnoise_f; //C noise in frequency space
fftwf_complex *Dnoise_f; //D noise in frequency space
fftwf_complex *channel1_f; //channel 1 record in frequency space
fftwf_complex *channel2_f; //channel 2 record in frequency space
fftwf_complex *signal_f; //signal without noise frequency space
float *channel1; //channel 1 record in time space
float *channel2; //channel 2 record in time space
float *signal_t; //signal without noise in time space
fftwf_plan channel1_plan;
fftwf_plan channel2_plan;
fftwf_plan signal_plan;
unsigned fft_flags=FFTW_ESTIMATE;

float total_mixing_frequency=25e9; //sum of local oscillator freqencies, in Hz
float sampling_rate=200e6; //digitizer sampling rate in Hz
float frequency_bin_width=sampling_rate/((float)record_size);

int nrecords=2; //how many records to generate

string eggname="test.egg"; //output egg file name
string mcinfoname="test.mcinfo"; //output mcinfo name

class ChirpEvent {
public:
    double start_energy; //in eV
    double duration; //in seconds
    double start_time; //in seconds
};

ChirpEvent *events;
int nevents=0;
double on_time=0; //progress in generation in seconds

double record_time=((double)record_size)/sampling_rate;

float getGaussianRand(float mean,float sigma);
float getUniformRand(float range);
void generate_record(unsigned char *dataptr);
void make_json_string_entry(yajl_gen gen,const char *key,const char *value);
void make_json_numeric_entry(yajl_gen gen,const char *key,double value);
void make_json_integer_entry(yajl_gen gen,const char *key,int value);
void yajl_gen_number_forreals(yajl_gen gen,double num);
void yajl_gen_integer_forreals(yajl_gen gen,int num);

int main(int argc,char *argv[])
{
    //allocate working space 
    Anoise_f=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*fft_size);
    Bnoise_f=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*fft_size);
    Cnoise_f=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*fft_size);
    Dnoise_f=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*fft_size);
    channel1_f=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*fft_size);
    channel2_f=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*fft_size);
    signal_f=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*fft_size);
    channel1=(float*)fftwf_malloc(sizeof(float)*record_size);
    channel2=(float*)fftwf_malloc(sizeof(float)*record_size);
    signal_t=(float*)fftwf_malloc(sizeof(float)*record_size);
    channel1_plan=fftwf_plan_dft_c2r_1d(record_size,channel1_f,channel1,fft_flags);
    channel2_plan=fftwf_plan_dft_c2r_1d(record_size,channel2_f,channel2,fft_flags);
    signal_plan=fftwf_plan_dft_r2c_1d(record_size,signal_t,signal_f,fft_flags);
    //Prep egg file
    Monarch *egg=Monarch::OpenForWriting(eggname.c_str());
    MonarchHeader *header=egg->GetHeader();
    header->SetAcqRate(sampling_rate/1e6);
    header->SetRecordSize(record_size);
    header->SetAcqTime((float)(record_size/sampling_rate)*nrecords);
    header->SetAcqMode(sTwoChannel);
    if(!egg->WriteHeader()) {
	cerr << "failed to write header" << endl;
    }
    //prep mcinfo file
    yajl_gen mcinfo_json=yajl_gen_alloc(NULL);
    yajl_gen_map_open(mcinfo_json);
    make_json_string_entry(mcinfo_json,"mc_id",mcname);
    make_json_numeric_entry(mcinfo_json,"system_temperature",20.4);
    make_json_string_entry(mcinfo_json,"noise_model","full_receiver_model");
    make_json_string_entry(mcinfo_json,"receiver_assumptions","TODO");

    //figure out events
    nevents=0;
//    events=new ChirpEvent[nevents];
    for(int i=0;i<nevents;i++) {
	events[i].start_time=0;
	events[i].duration=100;
	events[i].start_energy=emass*(BField*ecyclo/((total_mixing_frequency+50e6)*2.0*M_PI)-1);
    }

    //record event info
    yajl_gen_string(mcinfo_json,(unsigned char*)("events"),strlen("events"));
    yajl_gen_map_open(mcinfo_json);
    for(int i=0;i<nevents;i++) {
	//record support
	int startrecord=(int)(events[i].start_time/record_time);
	int startsamp=(int)((events[i].start_time-((double)startrecord)*record_time)*sampling_rate);
	int endrecord=(int)((events[i].start_time+events[i].duration)/record_time);
	int endsamp=(int)((events[i].start_time+events[i].duration-((double)endrecord)*record_time)*sampling_rate);
	if(endrecord>nrecords) {
	    endrecord=nrecords-1;
	    endsamp=record_size-1;
	}
	yajl_gen_string(mcinfo_json,(unsigned char*)("support"),strlen("support"));
	yajl_gen_array_open(mcinfo_json);
	yajl_gen_integer(mcinfo_json,startrecord);
	yajl_gen_integer(mcinfo_json,startsamp);
	yajl_gen_integer(mcinfo_json,endrecord);
	yajl_gen_integer(mcinfo_json,endsamp);
	yajl_gen_array_close(mcinfo_json);
	//record start energy
	make_json_numeric_entry(mcinfo_json,"start_energy",events[i].start_energy);
    }
    yajl_gen_map_close(mcinfo_json);

    //generate monte carlo
    for(int onrecord=0;onrecord<nrecords;onrecord++) {
	MonarchRecord *r=egg->GetRecordInterleaved();
	generate_record(r->fDataPtr);
	if(!egg->WriteRecord()) {
	    cerr << "failed to write record" << endl;
    	}
    }
    //save egg
    egg->Close();
    //save mcinfo
    yajl_gen_map_close(mcinfo_json);
    const unsigned char *buf;
    size_t len;
    yajl_gen_get_buf(mcinfo_json,&buf,&len);
    FILE *mcinfo_file=fopen(mcinfoname.c_str(),"w");
    fwrite(buf,1,len,mcinfo_file);
    fclose(mcinfo_file);
    yajl_gen_clear(mcinfo_json);
    yajl_gen_free(mcinfo_json);
    //free memory
    fftwf_free(Anoise_f);
    fftwf_free(Bnoise_f);
    fftwf_free(Cnoise_f);
    fftwf_free(Dnoise_f);
    fftwf_free(channel1_f);
    fftwf_free(channel2_f);
    fftwf_free(signal_f);
    fftwf_free(channel1);
    fftwf_free(channel2);
    fftwf_free(signal_t);
    fftwf_destroy_plan(channel1_plan);
    fftwf_destroy_plan(channel2_plan);
    fftwf_destroy_plan(signal_plan);
}

//assumes an interleaved data pointer
void generate_record(unsigned char *dataptr) 
{
    float R=50; //cable impedence, ohms
    float sigma_a=sqrt(R*kB*T_A/2.0);
    float sigma_b=sqrt(R*kB*T_B/2.0);
    float sigma_c=sqrt(R*kB*T_C/2.0);
    float sigma_d=sqrt(R*kB*T_D/2.0);


    //generate frequency space noise 
    for(int k=0;k<fft_size;k++) {
	Anoise_f[k]=getGaussianRand(0,sigma_a)+I*getGaussianRand(0,sigma_a);
	Bnoise_f[k]=getGaussianRand(0,sigma_b)+I*getGaussianRand(0,sigma_b);
	Cnoise_f[k]=getGaussianRand(0,sigma_c)+I*getGaussianRand(0,sigma_c);
	Dnoise_f[k]=getGaussianRand(0,sigma_d)+I*getGaussianRand(0,sigma_d);
    }

    //Add A+B+Ce^ilwg and (D+C+Be^ilwg)e^idl
    for(int i=0;i<fft_size;i++) {
	float f=total_mixing_frequency+frequency_bin_width*((float)i);
	float wg_phase=2*M_PI*f*lwg/c_wg;
	complex float wg_e=cexp(I*wg_phase);
	channel1_f[i]=Anoise_f[i]+Bnoise_f[i]+Cnoise_f[i]*wg_e;
	channel2_f[i]=Dnoise_f[i]+Cnoise_f[i]+Bnoise_f[i]*wg_e;
    }

    //inject signal
    bool signal_present=false;
    for(int i=0;i<record_size;i++) {
	signal_t[i]=0;
    }
    for(int k=0;k<nevents;k++) {
	if(events[k].start_time<on_time||(events[k].start_time+events[k].duration)>(on_time+record_time)) { //this event is in my range
	    signal_present=true;
	    double frequency=ecyclo*emass/(emass+events[k].start_energy);
	    cerr << "frequency is " << frequency << endl;
	    //TODO include energy loss here, putting 1e-15 watts for starters
	    double chirpiness=frequency*6250/emass;
	    chirpiness=0;
	    double vmax=sqrt(1e-15*R);
	    long start_sample=(events[k].start_time-on_time)*sampling_rate;
	    if(start_sample<0) start_sample=0;
	    long stop_sample=(events[k].start_time+events[k].duration-on_time)*sampling_rate;
	    if(stop_sample>=record_size) stop_sample=record_size-1;
//	    cerr << "start sample " << start_sample << " stop sample " << stop_sample << "record size " << record_size << endl;
	    double  fftnorm=1/((double)record_size);
	    for(int i=start_sample;i<=stop_sample;i++) {
		double t=on_time+((double)i)/sampling_rate;
		double dt=t-events[k].start_time;
		double phase=frequency*dt+chirpiness*dt*dt-total_mixing_frequency*dt;
		signal_t[i]+=fftnorm*vmax*cos(phase);
	    }
	}
    }
    if(signal_present) {
	fftwf_execute(signal_plan);
	for(int i=0;i<fft_size;i++) {
	    channel1_f[i]+=0.5*signal_f[i];
	    channel2_f[i]+=0.5*signal_f[i];
	}
    }
    
    //apply cable delays
    for(int i=0;i<fft_size;i++) {
	float f=total_mixing_frequency+frequency_bin_width*((float)i);
	float cable_phase=2*M_PI*f*dl/c_cable;
	complex float cable_e=cexp(I*cable_phase);
	channel2_f[i]=channel2_f[i]*cable_e;
    }

    //TODO apply receiver chain transfer function
    for(int i=0;i<fft_size;i++) {
	channel1_f[i]*=1e11;
	channel2_f[i]*=1e11;
    }

    //fft to construct the time series again
    fftwf_execute(channel1_plan);
    fftwf_execute(channel2_plan);

    //discretize
    float fft_scale=1/((float)record_size);
    for(int i=0;i<record_size;i++) {
	float x1=256*(fft_scale*channel1[i]/digitizer_fullscale+digitizer_fullscale);
	float x2=256*(fft_scale*channel2[i]/digitizer_fullscale+digitizer_fullscale);
	//clipping
	if(x1>255) x1=255;
	if(x2>255) x2=255;
	if(x1<0) x1=0;
	if(x2<0) x2=0;
	dataptr[2*i]=(unsigned char)x1;
	dataptr[2*i+1]=(unsigned char)x2;
    }
    on_time+=record_time;
}

float getUniformRand(float range)
{
    return range*((float)rand())/((float)RAND_MAX);
}

bool has_spare_grand=false;
float spare_grand;

float getGaussianRand(float mean,float sigma)
{
    //Polar form of the Box-Muller transformation
    //taken from http://www.taygeta.com/random/gaussian.html
    if(has_spare_grand) {has_spare_grand=false; return spare_grand*sigma+mean;}
   float x1, x2, w, y1;
   do {
             x1 = 2.0 * getUniformRand(1.0) - 1.0;
             x2 = 2.0 * getUniformRand(1.0) - 1.0;
             w = x1 * x1 + x2 * x2;
   } while ( w >= 1.0 );
   w = sqrt( (-2.0 * log( w ) ) / w );
   y1 = x1 * w;
    spare_grand= x2 * w;
    has_spare_grand=true;
   return y1*sigma+mean;
}


void make_json_string_entry(yajl_gen gen,const char *key,const char *value)
{
    yajl_gen_string(gen,(unsigned char*)key,strlen(key));
    yajl_gen_string(gen,(unsigned char*)value,strlen(value));
}

void make_json_numeric_entry(yajl_gen gen,const char *key,double value)
{
    yajl_gen_string(gen,(unsigned char*)key,strlen(key));
    yajl_gen_number_forreals(gen,value);
}

void make_json_integer_entry(yajl_gen gen,const char *key,int value)
{
    yajl_gen_string(gen,(unsigned char*)key,strlen(key));
    yajl_gen_integer(gen,value);
}


void yajl_gen_number_forreals(yajl_gen gen,double num)
{
    char buffer[256];
    sprintf(buffer,"%f",num);
    yajl_gen_number(gen,buffer,strlen(buffer));
}
    
