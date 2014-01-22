#include <complex.h>
#include "Monarch.hpp"
#include <math.h>
#include <fftw3.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include "ReceiverTransferFunctions.hh"
#include <yajl/yajl_gen.h>
#include <yajl/yajl_tree.h>
using namespace std;

//program constants
char mcname[]="locust_mc";

//physical constants
double kB=1.3806e-23;  // Boltzman constant W/K
double c=3e10; //speed of light in cm/s
double ecyclo=1.758820e11; //electron cyclotron frequency in rad/s
double emass=510998.9; //electron mass in eV
double JoulesToEv=1/1.6e-19; // eV/J

//     for setup DOUBLEAMP
// Waveguide is set up like this
// --<A|---<B>-----<C>----|D>-----
// ------------lwg--------------dl
// Where A,B,C,D are noise sources
// lwg and dl are the waveguide length and extra cable length
//
//     for setup SINGLEAMP
// Waveguide is set up like this
// --<A|---<B>----lwg-----trapregion----dl------REFLECTION
// where A is the noise of the amp that is not transmitted both ways
// B is the noise of the amp that is transmitted both directions
// lwg is the length of waveguide from the amplifier to the trapping region
// dl is the length of waveguide from the trap region to the short at the bottom
// there is a 100% reflective short at the bottom
// amplifier impedance mismatch is neglected

double c_wg=c*sqrt(1-pow(c/(2.0*0.42*2.54*26e9),2.0)); //speed of light in waveguide in cm/s
double c_cable=0.6*c; //speed of light in cables in cm/s (this is a guess)
double lwg=20; //waveguide length in cm (TODO measure this)
double dl=5; //extra length on channel 2 cable in cm
int record_size=4194304; //number of samples in record time series
int fft_size=record_size/2+1;
double digitizer_fullscale=0.5; //in volts
double amp_temp=30; //in Kelvin
double uncor_temp=30; //in Kelvin
double T_A=uncor_temp; //temperature of noise source A
double T_B=amp_temp; //temperature of noise source B
double T_C=amp_temp; //temperature of noise source C
double T_D=uncor_temp; //temperature of noise source D
double BField=0.90; //in tesla
fftw_complex *Anoise_f; //A noise in frequency space
fftw_complex *Bnoise_f; //B noise in frequency space
fftw_complex *Cnoise_f; //C noise in frequency space
fftw_complex *Dnoise_f; //D noise in frequency space
fftw_complex *channel1_f; //channel 1 record in frequency space
fftw_complex *channel2_f; //channel 2 record in frequency space
fftw_complex *signal_f; //signal without noise frequency space
double *channel1; //channel 1 record in time space
double *channel2; //channel 2 record in time space
double *signal_t; //signal without noise in time space
fftw_plan channel1_plan; //for fftw
fftw_plan channel2_plan; //for fftw
fftw_plan signal_plan; //for fftw
unsigned fft_flags=FFTW_ESTIMATE;

double hf_mixing_frequency=24.2e9; //high frequency oscillator in Hz
double lf_mixing_frequency=500e6;  //low frequency oscillator in Hz
double total_mixing_frequency=hf_mixing_frequency+lf_mixing_frequency; //sum of local oscillator freqencies, in Hz
double sampling_rate=200e6; //digitizer sampling rate in Hz
double frequency_bin_width=sampling_rate/((double)record_size); //how big one bin is, in Hz
string waveguide_setup="DOUBLEAMP";

int nrecords=10; //how many records to generate
double expected_event_rate=100;  //expected event rate in hz

string eggname="test.egg"; //output egg file name
string mcinfoname="test.mcinfo"; //output mcinfo name

//transfer function

string transfer_function_fname="receiver_transfer_functions.json";
ReceiverTransferFunctions tfuncs;
TransferFunction hf_transferfunction_ch1; //high frequency transfer function interpolated into current band
TransferFunction hf_transferfunction_ch2; //high frequency transfer function interpolated into current band
/*
double *hf_frequencies; //frequencies used for transfer function
double *hf_full_transferfunction_ch1; //transfer function NOT IN DB
double *hf_full_transferfunction_ch2;
int hf_arraycount; //size of big transfer function arrays
double *hf_transferfunction_ch1; //high frequency transfer function interpolated into current band
double *hf_transferfunction_ch2; //high frequency transfer function interpolated into current band
*/

class ChirpEvent {
public:
    double duration; //in seconds
    double start_time; //in seconds
    double start_frequency; //in Hz
    double dfdt; //in Hz/s
    double power; //power radiated into chirp in Watts
    //not actually used, but relevant to calculating above parameters
   // double start_energy; //in eV
};

ChirpEvent *events; //array of events to add
int nevents=10; //how many events are there
double on_time=0; //progress in generation in seconds

double record_time=((double)record_size)/sampling_rate; //length of one record in seconds

double getGaussianRand(double mean,double sigma);
double getUniformRand(double range);
double getExponentialRand(double scale);
void generate_record_singleamp(unsigned char *dataptr);
void generate_record_doubleamp(unsigned char *dataptr);
void make_json_string_entry(yajl_gen gen,const char *key,const char *value);
void make_json_numeric_entry(yajl_gen gen,const char *key,double value);
void make_json_integer_entry(yajl_gen gen,const char *key,int value);
void yajl_gen_number_forreals(yajl_gen gen,double num);
void yajl_gen_integer_forreals(yajl_gen gen,int num);
int load_transfer_functions(const char *fname);
double interpolate_transferfunction(double *freqs,double *vals,int len,double f);

double getJsonDouble(yajl_val node,const char *key);
string getJsonString(yajl_val node,const char *key);
int load_config_file(const char *fname);

int main(int argc,char *argv[])
{
    if(argc<2) {
	printf("usage: generate_data configfile\n");
	return -1;
    }
    //load the config file
    if(load_config_file(argv[1])) {
	printf("error reading config file");
	return -1;
    }
    //allocate working space 
    Anoise_f=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*fft_size);
    Bnoise_f=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*fft_size);
    Cnoise_f=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*fft_size);
    Dnoise_f=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*fft_size);
    channel1_f=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*fft_size);
    channel2_f=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*fft_size);
    signal_f=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*fft_size);
//    hf_transferfunction_ch1=new double[fft_size];
//    hf_transferfunction_ch2=new double[fft_size];
    channel1=(double*)fftw_malloc(sizeof(double)*record_size);
    channel2=(double*)fftw_malloc(sizeof(double)*record_size);
    signal_t=(double*)fftw_malloc(sizeof(double)*record_size);
    channel1_plan=fftw_plan_dft_c2r_1d(record_size,channel1_f,channel1,fft_flags);
    channel2_plan=fftw_plan_dft_c2r_1d(record_size,channel2_f,channel2,fft_flags);
    signal_plan=fftw_plan_dft_r2c_1d(record_size,signal_t,signal_f,fft_flags);
    //Load Receiver Transfer Function
    load_transfer_functions(transfer_function_fname.c_str());
    //
    //Prep egg file
    Monarch *egg=Monarch::OpenForWriting(eggname.c_str());
    MonarchHeader *header=egg->GetHeader();
    header->SetAcquisitionRate(sampling_rate/1e6);
    header->SetRecordSize(record_size);
    header->SetRunDuration(uint(nrecords*record_time));
//    header->SetAcquisitionTime((double)(record_size/sampling_rate)*nrecords);
    if(waveguide_setup=="DOUBLEAMP") {
        header->SetAcquisitionMode(sTwoChannel);
    } else if(waveguide_setup=="SINGLEAMP") {
        header->SetAcquisitionMode(sOneChannel);
    } else {
        cerr << "invalid waveguide setup: " << waveguide_setup << " should be DOUBLEAMP or SINGLEAMP" << endl;
        return -1;
    }
    egg->WriteHeader();
    //prep mcinfo file
    yajl_gen mcinfo_json=yajl_gen_alloc(NULL);
    yajl_gen_config(mcinfo_json, yajl_gen_beautify, 1);
    yajl_gen_map_open(mcinfo_json);
    make_json_string_entry(mcinfo_json,"mc_id",mcname);
    make_json_integer_entry(mcinfo_json,"record_size",header->GetRecordSize());
    //make_json_integer_entry(mcinfo_json,"records_simulated",);
    make_json_numeric_entry(mcinfo_json,"acquisition_rate",header->GetAcquisitionRate());
    make_json_numeric_entry(mcinfo_json,"system_temperature",90);
    make_json_string_entry(mcinfo_json,"noise_model","full_receiver_model_minus_lf_mixing");
    make_json_string_entry(mcinfo_json,"receiver_assumptions","amplifiers correlated");
    make_json_string_entry(mcinfo_json,"egg_name",eggname.c_str());
    make_json_numeric_entry(mcinfo_json,"lf_mixing_frequency",lf_mixing_frequency);
    make_json_numeric_entry(mcinfo_json,"hf_mixing_frequency",hf_mixing_frequency);

    //place some events
    /*  //this is depricated.  Events are provided in config file now
    nevents=5;
    events=new ChirpEvent[nevents];
    for(int i=0;i<nevents;i++) {
	events[i].start_time=getUniformRand(record_time*((double)nrecords));
	events[i].duration=getExponentialRand(200e-6);
	events[i].start_energy=emass*(BField*ecyclo/((total_mixing_frequency+50e6)*2.0*M_PI)-1);
	events[i].start_frequency=ecyclo*emass/(emass+events[i].start_energy);
	events[i].power=1e-15;
	events[i].dfdt=events[i].start_frequency*events[i].power*JoulesToEv/emass;
    }
    */

    //record event info in json file
    yajl_gen_string(mcinfo_json,(unsigned char*)("events"),strlen("events"));
    yajl_gen_array_open(mcinfo_json);
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
    	yajl_gen_map_open(mcinfo_json);
    	yajl_gen_string(mcinfo_json,(unsigned char*)("support"),strlen("support"));
    	yajl_gen_array_open(mcinfo_json);
    	yajl_gen_integer(mcinfo_json,startrecord);
    	yajl_gen_integer(mcinfo_json,startsamp);
    	yajl_gen_integer(mcinfo_json,endrecord);
    	yajl_gen_integer(mcinfo_json,endsamp);
    	yajl_gen_array_close(mcinfo_json);
    	//record start energy
    	//make_json_numeric_entry(mcinfo_json,"start_energy",events[i].start_energy);
    	make_json_numeric_entry(mcinfo_json,"start_frequency",events[i].start_frequency);
    	make_json_numeric_entry(mcinfo_json,"power",events[i].power);
    	make_json_numeric_entry(mcinfo_json,"dfdt",events[i].dfdt);
    	make_json_numeric_entry(mcinfo_json,"duration",events[i].duration);
    	yajl_gen_map_close(mcinfo_json);
    }
    yajl_gen_array_close(mcinfo_json);

    //generate monte carlo
    for(int onrecord=0;onrecord<nrecords;onrecord++) {
        if(waveguide_setup=="DOUBLEAMP") {
    	    MonarchRecord *r=egg->GetRecordInterleaved();
        	generate_record_doubleamp(r->fData);
        } else {
    	    MonarchRecord *r=egg->GetRecordSeparateOne();
        	generate_record_singleamp(r->fData);
        }
    	if(!egg->WriteRecord()) {
    	    cerr << "failed to write record" << endl;
        	}
    }
    make_json_integer_entry(mcinfo_json,"records_simulated",nrecords);
    
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
    fftw_free(Anoise_f);
    fftw_free(Bnoise_f);
    fftw_free(Cnoise_f);
    fftw_free(Dnoise_f);
    fftw_free(channel1_f);
    fftw_free(channel2_f);
    fftw_free(signal_f);
    fftw_free(channel1);
    fftw_free(channel2);
    fftw_free(signal_t);
    fftw_destroy_plan(channel1_plan);
    fftw_destroy_plan(channel2_plan);
    fftw_destroy_plan(signal_plan);
}

//assumes an interleaved data pointer
void generate_record_singleamp(unsigned char *dataptr) 
{
    double R=50; //cable impedence, ohms (so I can turn power into volts)
    double sigma_a=sqrt(R*kB*T_A/2.0);
    double sigma_b=sqrt(R*kB*T_B/2.0);

    double total_length=2.0*(lwg+dl); //total length from amp to bottom and back

    //generate frequency space noise 
    for(int k=0;k<fft_size;k++) {
	Anoise_f[k]=getGaussianRand(0,sigma_a)+I*getGaussianRand(0,sigma_a);
	Bnoise_f[k]=getGaussianRand(0,sigma_b)+I*getGaussianRand(0,sigma_b);
    }

    //combine noise sources
    //Add A+B+Ce^ilwg and (D+C+Be^ilwg)e^idl
    for(int i=0;i<fft_size;i++) {
	double f=total_mixing_frequency+frequency_bin_width*((double)i);
	double wg_phase=2*M_PI*f*total_length/c_wg;
	complex double wg_e=cexp(I*wg_phase);
	channel1_f[i]=Anoise_f[i]+Bnoise_f[i]+Bnoise_f[i]*wg_e;
    }

    //--inject signal--
    //first make signal in time space
    bool signal_present=false;
    for(int i=0;i<record_size;i++) {
	signal_t[i]=0;
    }
    printf("on time %g to %g\n",on_time,on_time+record_time);
    for(int k=0;k<nevents;k++) {
//	if(events[k].start_time<on_time||(events[k].start_time+events[k].duration)>(on_time+record_time)) { //this event is in my range
	if((events[k].start_time<(on_time+record_time))&&((events[k].start_time+events[k].duration)>(on_time))) { //this event is in my range
        printf("creating event\n");
	    signal_present=true;
	    double vmax=sqrt(events[k].power*R);
	    long start_sample=(events[k].start_time-on_time)*sampling_rate;
	    if(start_sample<0) start_sample=0;
	    long stop_sample=(events[k].start_time+events[k].duration-on_time)*sampling_rate;
	    if(stop_sample>=record_size) stop_sample=record_size-1;
	    double  fftnorm=1/((double)record_size);
	    for(int i=start_sample;i<=stop_sample;i++) {
		double t=on_time+((double)i)/sampling_rate;
		double dt=t-events[k].start_time;
		double phase=(events[k].start_frequency-total_mixing_frequency)*dt+events[k].dfdt*dt*dt;
		signal_t[i]+=fftnorm*vmax*cos(2.0*M_PI*phase);
	    }
	}
    }
    //now turn signal into frequency space
    if(signal_present) {
	fftw_execute(signal_plan);
	for(int i=0;i<fft_size;i++) {
    	double f=total_mixing_frequency+frequency_bin_width*((double)i);
    	double wg_phase=2*M_PI*f*2.0*dl/c_wg;
	    complex double wg_e=cexp(I*wg_phase);
	    channel1_f[i]+=0.5*signal_f[i]+0.5*wg_e*signal_f[i];
//	    channel2_f[i]+=0.5*signal_f[i];
	}
    }

    //cable delays don't matter, ignore
    //
    //apply high and low frequency transfer function
    for(int i=0;i<fft_size;i++) {
	  channel1_f[i]*=hf_transferfunction_ch1.data[i];
	//channel2_f[i]*=hf_transferfunction_ch2[i];
    }

    //TODO apply low frequency receiver chain transfer function
    
    for(int i=0;i<fft_size;i++) {
	//this applies uniform gain
	//as a stopgap measure until real gain is known
	channel1_f[i]*=4e13;
//	channel2_f[i]*=1e7;
    }
   

    //fft to construct the time series again
    fftw_execute(channel1_plan);

    //discretize
    double fft_scale=1/((double)record_size);
    for(int i=0;i<record_size;i++) {
	double x1=256*(fft_scale*channel1[i]/digitizer_fullscale+digitizer_fullscale);
	double x2=256*(fft_scale*channel2[i]/digitizer_fullscale+digitizer_fullscale);
	//clipping to 8 bits
	if(x1>255) x1=255;
	if(x2>255) x2=255;
	if(x1<0) x1=0;
	if(x2<0) x2=0;
    //don't skip every other one, we're only saving one channel
	dataptr[i]=(unsigned char)x1;
//	dataptr[2*i+1]=(unsigned char)x2;
    }
    //update time in seconds
    on_time+=record_time;
}

void generate_record_doubleamp(unsigned char *dataptr) 
{
    double R=50; //cable impedence, ohms (so I can turn power into volts)
    double sigma_a=sqrt(R*kB*T_A/2.0);
    double sigma_b=sqrt(R*kB*T_B/2.0);
    double sigma_c=sqrt(R*kB*T_C/2.0);
    double sigma_d=sqrt(R*kB*T_D/2.0);
    cerr << "sigma_a " << sigma_a << endl;
    cerr << "sigma_b " << sigma_b << endl;
    cerr << "sigma_c " << sigma_c << endl;
    cerr << "sigma_d " << sigma_d << endl;


    //generate frequency space noise 
    for(int k=0;k<fft_size;k++) {
	Anoise_f[k]=getGaussianRand(0,sigma_a)+I*getGaussianRand(0,sigma_a);
	Bnoise_f[k]=getGaussianRand(0,sigma_b)+I*getGaussianRand(0,sigma_b);
	Cnoise_f[k]=getGaussianRand(0,sigma_c)+I*getGaussianRand(0,sigma_c);
	Dnoise_f[k]=getGaussianRand(0,sigma_d)+I*getGaussianRand(0,sigma_d);
    }

    //combine noise sources
    //Add A+B+Ce^ilwg and (D+C+Be^ilwg)e^idl
    for(int i=0;i<fft_size;i++) {
	double f=total_mixing_frequency+frequency_bin_width*((double)i);
	double wg_phase=2*M_PI*f*lwg/c_wg;
	complex double wg_e=cexp(I*wg_phase);
	channel1_f[i]=Anoise_f[i]+Bnoise_f[i]+Cnoise_f[i]*wg_e;
	channel2_f[i]=Dnoise_f[i]+Cnoise_f[i]+Bnoise_f[i]*wg_e;
    }

    //--inject signal--
    //first make signal in time space
    bool signal_present=false;
    for(int i=0;i<record_size;i++) {
	signal_t[i]=0;
    }
    printf("on time %g to %g\n",on_time,on_time+record_time);
    for(int k=0;k<nevents;k++) {
//	if(events[k].start_time<on_time||(events[k].start_time+events[k].duration)>(on_time+record_time)) { //this event is in my range
	if((events[k].start_time<(on_time+record_time))&&((events[k].start_time+events[k].duration)>(on_time))) { //this event is in my range
        printf("creating event\n");
	    signal_present=true;
	    double vmax=sqrt(events[k].power*R);
	    long start_sample=(events[k].start_time-on_time)*sampling_rate;
	    if(start_sample<0) start_sample=0;
	    long stop_sample=(events[k].start_time+events[k].duration-on_time)*sampling_rate;
	    if(stop_sample>=record_size) stop_sample=record_size-1;
	    double  fftnorm=1/((double)record_size);
	    for(int i=start_sample;i<=stop_sample;i++) {
		double t=on_time+((double)i)/sampling_rate;
		double dt=t-events[k].start_time;
		double phase=(events[k].start_frequency-total_mixing_frequency)*dt+events[k].dfdt*dt*dt;
		signal_t[i]+=fftnorm*vmax*cos(2.0*M_PI*phase);
	    }
	}
    }
    //now turn signal into frequency space
    if(signal_present) {
	fftw_execute(signal_plan);
	for(int i=0;i<fft_size;i++) {
	    channel1_f[i]+=0.5*signal_f[i];
	    channel2_f[i]+=0.5*signal_f[i];
	}
    }
    
    //apply cable delays
    for(int i=0;i<fft_size;i++) {
	double f=total_mixing_frequency+frequency_bin_width*((double)i);
	double cable_phase=2*M_PI*f*dl/c_cable;
	complex double cable_e=cexp(I*cable_phase);
	channel2_f[i]=channel2_f[i]*cable_e;
    }

    //apply high and low frequency transfer function
    /*
    for(int i=0;i<fft_size;i++) {
	channel1_f[i]*=hf_transferfunction_ch1.data[i];
	channel2_f[i]*=hf_transferfunction_ch2.data[i];
    }
    */

    //TODO apply low frequency receiver chain transfer function
    for(int i=0;i<fft_size;i++) {
	//this applies uniform 70 dB gain
	//as a stopgap measure until real gain is known
	channel1_f[i]*=1e7;
	channel2_f[i]*=1e7;
    }

    //fft to construct the time series again
    fftw_execute(channel1_plan);
    fftw_execute(channel2_plan);

    //discretize
    double fft_scale=1/((double)record_size);
    for(int i=0;i<record_size;i++) {
	double x1=256*(fft_scale*channel1[i]/digitizer_fullscale+digitizer_fullscale);
	double x2=256*(fft_scale*channel2[i]/digitizer_fullscale+digitizer_fullscale);
	//clipping to 8 bits
	if(x1>255) x1=255;
	if(x2>255) x2=255;
	if(x1<0) x1=0;
	if(x2<0) x2=0;
	dataptr[2*i]=(unsigned char)x1;
	dataptr[2*i+1]=(unsigned char)x2;
    }
    //update time in seconds
    on_time+=record_time;
}

double getUniformRand(double range)
{
    return range*((double)rand())/((double)RAND_MAX);
}

double getExponentialRand(double scale)
{
    return -scale*log(getUniformRand(1.0));
}

bool has_spare_grand=false;
double spare_grand;

double getGaussianRand(double mean,double sigma)
{
    //Polar form of the Box-Muller transformation
    //taken from http://www.taygeta.com/random/gaussian.html
    if(has_spare_grand) {has_spare_grand=false; return spare_grand*sigma+mean;}
   double x1, x2, w, y1;
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


//-----some json helper functions---
//
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
    sprintf(buffer,"%g",num);
    yajl_gen_number(gen,buffer,strlen(buffer));
}

int load_transfer_functions(const char *fname)
{
    tfuncs.load_from_json(string(fname));
    hf_transferfunction_ch1=tfuncs.getTransferFunction(lf_mixing_frequency/1e6,sampling_rate/2/1e6,record_size/2);
    hf_transferfunction_ch2=tfuncs.getTransferFunction(lf_mixing_frequency/1e6,sampling_rate/2/1e6,record_size/2);
    /*
    //I just copied this from the yajl example
    //it loads the transfer function in json form and then interpolates around
    //the desired bandwidth
    FILE *file=fopen(fname,"r");
    char errbuf[1024];
    static unsigned char fileData[65536];
    fileData[0] = errbuf[0] = 0;
    size_t rd = fread((void *) fileData, 1, sizeof(fileData) - 1, file);
    if (rd == 0 && !feof(stdin)) {
        fprintf(stderr, "error encountered on file read\n");
        return 1;
    } else if (rd >= sizeof(fileData) - 1) {
        fprintf(stderr, "config file too big\n");
        return 1;
    }
    yajl_val node;
    node = yajl_tree_parse((const char *) fileData, errbuf, sizeof(errbuf));
    if (node == NULL) {
        fprintf(stderr, "parse_error: ");
        if (strlen(errbuf)) fprintf(stderr, " %s", errbuf);
        else fprintf(stderr, "unknown error");
        fprintf(stderr, "\n");
        return 1;
    }   
    const char *freqpath[]={"hf_frequencies",(const char*)0};
    const char *hf1path[]={"hf_transfer_channel1",(const char*)0};
    const char *hf2path[]={"hf_transfer_channel2",(const char*)0};
    yajl_val freqs=yajl_tree_get(node,freqpath,yajl_t_array);
    yajl_val hf1=yajl_tree_get(node,hf1path,yajl_t_array);
    yajl_val hf2=yajl_tree_get(node,hf2path,yajl_t_array);
    hf_frequencies=new double[freqs->u.array.len];
    hf_full_transferfunction_ch1=new double[freqs->u.array.len];
    hf_full_transferfunction_ch2=new double[freqs->u.array.len];
    hf_arraycount=freqs->u.array.len;
    for(size_t i=0;i<freqs->u.array.len;i++) {
	hf_frequencies[i]=freqs->u.array.values[i]->u.number.d;
	hf_full_transferfunction_ch1[i]=pow10(0.1*hf1->u.array.values[i]->u.number.d);
	hf_full_transferfunction_ch2[i]=pow10(0.1*hf2->u.array.values[i]->u.number.d);
    }
    yajl_tree_free(node);
    //now create the transfer function for my bandwidth
    for(int i=0;i<fft_size;i++) {
	double f=total_mixing_frequency+frequency_bin_width*((double)i);
	hf_transferfunction_ch1[i]=interpolate_transferfunction(hf_frequencies,hf_full_transferfunction_ch1,hf_arraycount,f);
	hf_transferfunction_ch2[i]=interpolate_transferfunction(hf_frequencies,hf_full_transferfunction_ch2,hf_arraycount,f);
    }
    */
    return 0;
}

//linear interpolation between frequency points
double interpolate_transferfunction(double *freqs,double *vals,int len,double f)
{
    if(f<=freqs[0]) return vals[0];
    if(f>=freqs[len-1]) return vals[len-1];
    for(int i=0;i<len;i++) {
	if( (f>=freqs[i])&&(f<=freqs[i+1])) {
	    double x1=freqs[i];
	    double x2=freqs[i+1];
	    double y1=vals[i];
	    double y2=vals[i+1];
	    return y1+(y2-y1)*((f-x1)/(x2-x1));
	}
    }
    cerr << "interpolation failed" << endl;
    return nanf("bad programmer");
}

int load_config_file(const char *fname) {
    //I just copied this from the yajl example
    //it loads the transfer function in json form and then interpolates around
    //the desired bandwidth
    FILE *file=fopen(fname,"r");
    char errbuf[1024];
    static unsigned char fileData[65536*16];
    fileData[0] = errbuf[0] = 0;
    size_t rd = fread((void *) fileData, 1, sizeof(fileData) - 1, file);
    if (rd == 0 && !feof(stdin)) {
        fprintf(stderr, "error encountered on file read\n");
        return 1;
    } else if (rd >= sizeof(fileData) - 1) {
        fprintf(stderr, "config file too big\n");
        return 1;
    }
    yajl_val node;
    node = yajl_tree_parse((const char *) fileData, errbuf, sizeof(errbuf));
    if (node == NULL) {
        fprintf(stderr, "parse_error: ");
        if (strlen(errbuf)) fprintf(stderr, " %s", errbuf);
        else fprintf(stderr, "unknown error");
        fprintf(stderr, "\n");
        return 1;
        }
    waveguide_setup=getJsonString(node,"waveguide_setup");
    transfer_function_fname=getJsonString(node,"transfer_function_filename");
    if(waveguide_setup=="DOUBLEAMP") {
        T_A=getJsonDouble(node,"receiver1_noise_temperature");
        T_B=getJsonDouble(node,"amp1_noise_temperature");
        T_C=getJsonDouble(node,"amp2_noise_temperature");
        T_D=getJsonDouble(node,"receiver2_noise_temperature");
        lwg=getJsonDouble(node,"waveguide_length");
        dl=getJsonDouble(node,"phase_delay_length");
    } else if (waveguide_setup=="SINGLEAMP") {
        T_A=getJsonDouble(node,"receiver1_noise_temperature");
        T_B=getJsonDouble(node,"amp1_noise_temperature");
        lwg=getJsonDouble(node,"waveguide_length");
        dl=getJsonDouble(node,"distance_to_short");
    }
    BField=getJsonDouble(node,"BField");
    hf_mixing_frequency=getJsonDouble(node,"hf_mixing_frequency");
    lf_mixing_frequency=getJsonDouble(node,"lf_mixing_frequency");
    total_mixing_frequency=hf_mixing_frequency+lf_mixing_frequency;
    double total_duration=getJsonDouble(node,"datafile_duration");
    nrecords=(int)(ceil(total_duration/record_time));
    eggname=getJsonString(node,"egg_outfile_name");
    mcinfoname=getJsonString(node,"mcinfo_outfile_name");
    //load the events
    const char *eventpath[]={"events",(const char *)0};
    yajl_val yevents=yajl_tree_get(node,eventpath,yajl_t_array);
    nevents=yevents->u.array.len;
    events=new ChirpEvent[nevents];
    for(int i=0;i<nevents;i++) {
	events[i].start_time=getJsonDouble(yevents->u.array.values[i],"start_time");
	events[i].duration=getJsonDouble(yevents->u.array.values[i],"duration");
	events[i].start_frequency=getJsonDouble(yevents->u.array.values[i],"start_frequency");
	events[i].dfdt=getJsonDouble(yevents->u.array.values[i],"dfdt");
	events[i].power=getJsonDouble(yevents->u.array.values[i],"power");
    }

    return 0;
}

double getJsonDouble(yajl_val node,const char *key)
{
    //const char *freqpath[]={"hf_frequencies",(const char*)0};
    const char *path[]={key,(const char*)0};
    yajl_val val=yajl_tree_get(node,path,yajl_t_number);
    if(val==NULL) {
	cerr << "error getting json double " << key << endl;
	return 0;
    }
    return val->u.number.d;
}

string getJsonString(yajl_val node,const char *key)
{
    //const char *freqpath[]={"hf_frequencies",(const char*)0};
    const char *path[]={key,(const char*)0};
    yajl_val val=yajl_tree_get(node,path,yajl_t_string);
    if(val==NULL) {
	cerr << "error getting json string " << key << endl;
	return "";
    }
    return string(val->u.string);
}
