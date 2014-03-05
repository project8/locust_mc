#include "ReceiverTransferFunctions.hh"
#include <math.h>

void tyajl_gen_string(yajl_gen node,string s);
double tgetJsonDouble(yajl_val node,const char *key);
void tyajl_gen_numberforreals(yajl_gen gen,double num);
void tmake_jsan_numeric_entry(yajl_gen gen,const char *key,double value);

TransferFunction::TransferFunction()
{
    data=NULL;
    size=0;
}

TransferFunction::TransferFunction(double fstart,double fspan,unsigned int sz)
{
    data=NULL;
    resize(sz);
    frequency_start=fstart;
    frequency_span=fspan;
}

TransferFunction::TransferFunction(const TransferFunction &f)
{
    data=NULL;
    (*this)=f;
}

TransferFunction::~TransferFunction()
{
    delete data;
}

void TransferFunction::resize(unsigned sz)
{
    delete data;
    size=sz;
    data=new float[size];
}
    
TransferFunction &TransferFunction::operator=(const TransferFunction &f)
{
    resize(f.size);
    frequency_start=f.frequency_start;
    frequency_span=f.frequency_span;
    memcpy(data,f.data,sizeof(float)*size);
    return *this;
}

TransferFunction TransferFunction::interpolateTo(double fstart,double fspan,unsigned int sz)
{
    TransferFunction ret(fstart,fspan,sz);
    for(unsigned int i=0;i<sz;i++) {
        double f=fstart+fspan*((double)i)/((double)sz);
        ret.data[i]=interpolatePoint(f);
    }
    return ret;
}
    
float TransferFunction::interpolatePoint(double freq)
{
    int bottom=floor(((double)size)*(freq-frequency_start)/frequency_span);
    int top=bottom+1;
    if(bottom>=(int)size) return data[size-1];
    if(top<=0) return data[0];
    double bottompos=frequency_start+frequency_span*((double)bottom)/((double)size);
    double toppos=frequency_start+frequency_span*((double)top)/((double)size);
    return data[bottom]+(data[top]-data[bottom])*(freq-bottompos)/(toppos-bottompos);
}
    
void TransferFunction::loadFromFile(string fname)
{
    ifstream fin(fname.c_str());
    if(!fin.good()) {
        cerr << "unable to open " << fname << endl;
        return;
    }
    //count lines
    string line;
    unsigned int nlines=0;
    while(getline(fin,line)) {
        nlines++;
    }
    size=nlines;
    data=new float[size];
    fin.clear(); //clear end of fileflag
    fin.seekg(ios_base::beg);
    //get data
    float pow;
    float freq;
    unsigned int online=0;
    while(getline(fin,line)) {
        stringstream ss(line);
        ss >> freq >> pow;
        data[online]=pow;
        online++;
    }
    frequency_span=freq;


    /*
    cout << "on file " << fname << endl;
    FILE *file=fopen(fname.c_str(),"r");
    char errbuf[1024];

    static unsigned char fileData[65536*16];
    fileData[0] = errbuf[0] = 0;
    size_t rd = fread((void *) fileData, 1, sizeof(fileData) - 1, file);
    if(rd==0) {
        cerr << "error encountered reading file " << fname << endl;
        return;
    } else if(rd>=sizeof(fileData)-1) {
        cerr << "file " << fname << " too big!" << endl;
        return;
    }
    cout << " first character is *" << fileData[0] << "*" << endl;

    yajl_val node;
    node = yajl_tree_parse((const char *) fileData, errbuf, sizeof(errbuf));
if (node == NULL) {
        fprintf(stderr, "parse_error: ");
        if (strlen(errbuf)) fprintf(stderr, " %s", errbuf);
        else fprintf(stderr, "unknown error");
        fprintf(stderr, "\n");
        return ;
    }


    frequency_span=tgetJsonDouble(node,"sampling_rate")/2.0;
    const char *addr[]={"data",(const char*)0};
    yajl_val ydat=yajl_tree_get(node,addr,yajl_t_array);
    size=ydat->u.array.len;
    data=new float[ydat->u.array.len];
    for(unsigned int i=0;i<ydat->u.array.len;i++) {
        data[i]=ydat->u.array.values[i]->u.number.d;
    }
    yajl_tree_free(node);
    fclose(file);
    */
}
    
void TransferFunction::saveToJson(yajl_gen node)
{
    yajl_gen_map_open(node);
    tmake_jsan_numeric_entry(node,"frequency_start",frequency_start);
    tmake_jsan_numeric_entry(node,"frequency_span",frequency_span);
    tyajl_gen_string(node,"data");
    yajl_gen_array_open(node);
    for(unsigned int i=0;i<size;i++) {
        tyajl_gen_numberforreals(node,data[i]);
    }
    yajl_gen_array_close(node);
    yajl_gen_map_close(node);
}
    
void TransferFunction::loadFromJson(yajl_val node)
{
    delete data;
    const char *fstart_addr[]={"frequency_start",(const char*)0};
    const char *fspan_addr[]={"frequency_span",(const char*)0};
    const char *data_addr[]={"data",(const char*)0};
    yajl_val fstart_node=yajl_tree_get(node,fstart_addr,yajl_t_number);
    yajl_val fspan_node=yajl_tree_get(node,fspan_addr,yajl_t_number);
    yajl_val data_node=yajl_tree_get(node,data_addr,yajl_t_array);
    frequency_start=fstart_node->u.number.d;
    frequency_span=fspan_node->u.number.d;
    size=data_node->u.array.len;
    data=new float[size];
    for(size_t i=0;i<size;i++) {
        data[i]=data_node->u.array.values[i]->u.number.d;
    }
}
    
//--------------------ReceiverTransferFunctions----------

void ReceiverTransferFunctions::load_lf_from_filelist(string fname)
{
    ifstream fin(fname.c_str());
    string line;
    while(getline(fin,line)) {
        if(line.size()==0 ||line[0]=='#') continue;
        stringstream ss(line);
        float freq;
        string fname;
        ss >> freq >> fname;
        TransferFunction newfunc;
//        cout << "loading " <<  fname << endl;
        newfunc.loadFromFile(fname);
        newfunc.frequency_start=freq;
        low_frequency_stage[freq]=newfunc;
  //      cout << "size " << low_frequency_stage[freq].size << endl;
    }
  //  cout << "loaded " << low_frequency_stage.size() << " funcs" << endl;
    fin.close();
}
    
void ReceiverTransferFunctions::save_to_file(string fname)
{
    yajl_gen yout=yajl_gen_alloc(NULL);
    yajl_gen_config(yout,yajl_gen_beautify,1);
    yajl_gen_map_open(yout);
    tyajl_gen_string(yout,"hf_func");
    high_frequency_stage.saveToJson(yout);
    tyajl_gen_string(yout,"lf_funcs");
    yajl_gen_array_open(yout);
    for(map<float,TransferFunction>::iterator it=low_frequency_stage.begin();it!=low_frequency_stage.end();it++) {
        cout << "saving " << (*it).first << endl;
        (*it).second.saveToJson(yout);
    }
    yajl_gen_array_close(yout);
    yajl_gen_map_close(yout);
    size_t len;
    const unsigned char *buf;
    yajl_gen_get_buf(yout,&buf,&len);
    FILE *outfile=fopen(fname.c_str(),"w");
    fwrite(buf,1,len,outfile);
    fclose(outfile);
    yajl_gen_clear(yout);
    yajl_gen_free(yout);

}

double tgetJsonDouble(yajl_val node,const char *key)
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

void tmake_jsan_numeric_entry(yajl_gen gen,const char *key,double value)
{
    yajl_gen_string(gen,(unsigned char*)key,strlen(key));
    tyajl_gen_numberforreals(gen,value);
}

void tyajl_gen_numberforreals(yajl_gen gen,double num)
{
    char buffer[256];
    sprintf(buffer,"%g",num);
    yajl_gen_number(gen,buffer,strlen(buffer));
}


void tyajl_gen_string(yajl_gen node,string s)
{
    yajl_gen_string(node,(const unsigned char*)s.c_str(),s.size());
}
    
void ReceiverTransferFunctions::load_from_json(string fname)
{
    FILE *file=fopen(fname.c_str(),"r");
    char errbuf[1024];
    static unsigned char fileData[65536*16];
    fileData[0] = errbuf[0] = 0;
    size_t rd = fread((void *) fileData, 1, sizeof(fileData) - 1, file);
    if(rd==0) {
        cerr << "error encountered reading file " << fname << endl;
        return;
    } else if(rd>=sizeof(fileData)-1) {
        cerr << "file " << fname << " too big!" << endl;
        return;
    }
    yajl_val node;
    node = yajl_tree_parse((const char *) fileData, errbuf, sizeof(errbuf));
if (node == NULL) {
        fprintf(stderr, "parse_error: ");
        if (strlen(errbuf)) fprintf(stderr, " %s", errbuf);
        else fprintf(stderr, "unknown error");
        fprintf(stderr, "\n");
        return ;
    }
    const char *hfreq_addr[]={"hf_func",(const char*)0};
    const char *lfreqs_addr[]={"lf_funcs",(const char*)0};
    yajl_val hf_node=yajl_tree_get(node,hfreq_addr,yajl_t_object);
    yajl_val lfs_node=yajl_tree_get(node,lfreqs_addr,yajl_t_array);
    high_frequency_stage.loadFromJson(hf_node);
    for(size_t i=0;i<lfs_node->u.array.len;i++) {
        TransferFunction thefunc;
        thefunc.loadFromJson(lfs_node->u.array.values[i]);
        low_frequency_stage[thefunc.frequency_start]=thefunc;
    }
    yajl_tree_free(node);
    fclose(file);
}
    
TransferFunction ReceiverTransferFunctions::getTransferFunction(double lo_freq,double fspan,unsigned int sz)
{
    //get the transfer function from  high frequency stage
    TransferFunction hf=high_frequency_stage.interpolateTo(lo_freq+hf_lo_freq,fspan,sz);
    //find which lf transfer functions bracket my desired frequency
    map<float,TransferFunction>::iterator top=low_frequency_stage.upper_bound(lo_freq);
    TransferFunction lf;
    if(top==low_frequency_stage.end()) {
        map<float,TransferFunction>::reverse_iterator atend=(low_frequency_stage.rbegin());
        lf=(*atend).second.interpolateTo( (*atend).second.frequency_start,fspan,sz);
    } else {
        if(top==low_frequency_stage.begin()) {
            lf=(*top).second.interpolateTo((*top).second.frequency_start,fspan,sz);
        } else {
            map<float,TransferFunction>::iterator bottom=top;
            bottom--;
            TransferFunction lfbottom=(*bottom).second.interpolateTo((*bottom).second.frequency_start,fspan,sz);
            TransferFunction lftop=(*top).second.interpolateTo((*top).second.frequency_start,fspan,sz);
            double x1=(*bottom).second.frequency_start;
            double x2=(*top).second.frequency_start;
            double x=lo_freq;
            lf=TransferFunction(0,fspan,sz);
            for(unsigned int i=0;i<lf.size;i++) {
                float y1=lfbottom.data[i];
                float y2=lftop.data[i];
                lf.data[i]=y1+(y2-y1)*(x-x1)/(x2-x1);
            }
        }
    }
    TransferFunction ret(0,fspan,sz);
    for(unsigned int i=0;i<lf.size;i++) {
        ret.data[i]=hf.data[i]*lf.data[i];
    }
    return ret;
}
