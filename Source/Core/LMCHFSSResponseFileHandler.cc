/*
 * LMCHFSSFileHandler.cc
 *
 *  Created on: May 11, 2018
 *      Author: P. T. Surukuchi
 */

#include "LMCConst.hh"
#include "LMCException.hh"
#include "LMCHFSSResponseFileHandler.hh"

#include <stdio.h>
#include "logger.hh"

#include<unistd.h>

namespace locust
{
    LOGGER( lmclog, "HFSSResponseFileHandlerCore" );
    
    HFSSResponseFileHandlerCore::HFSSResponseFileHandlerCore():
    fTFNBins(1000),
    fFIRNBins(2000),
    fCropIndex(2100),
    fResolution(1e-12),
    fCharacteristicImpedance(1.0),
    fNSkips(1),
    fComplexFFT(),
    fHFSSFiletype(""),
    fIsFIRCreated(false),
    fWindowName("tukey"),
    fWindowParam(0.5),
    fPrintFIR ( false ),
    fConvertStoZ ( false ),
    fOutputPath( TOSTRING(PB_OUTPUT_DIR))
    {
    }
    
    HFSSResponseFileHandlerCore::~HFSSResponseFileHandlerCore()
    {
    }
    
    bool HFSSResponseFileHandlerCore::Configure(const scarab::param_node& aParam)
    {

        if( aParam.has( "print-fir-debug" ) )
    	{
    	    fPrintFIR = aParam["print-fir-debug"]().as_bool();
    	}

        if( aParam.has( "convert-sparams-to-z"))
        {
        	fConvertStoZ = aParam["convert-sparams-to-z"]().as_bool();
        }

        if( aParam.has( "crop-index"))
        {
        	fCropIndex = aParam["crop-index"]().as_int();
        }

        if( aParam.has( "characteristic-impedance"))
        {
            fCharacteristicImpedance = aParam["characteristic-impedance"]().as_double();
            LPROG( lmclog, "Characteristic impedance has been changed to " << fCharacteristicImpedance);
        }

    	if ( aParam.has( "output-path" ) )
    	{
    		fOutputPath = aParam["output-path"]().as_string();
    	}

        return true;
    }
    
    bool HFSSResponseFileHandlerCore::ends_with(const std::string &str, const std::string &suffix)
    {
        //copied from https://stackoverflow.com/a/20446239
        return str.size() >= suffix.size() &&
        str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
    }
    
    bool HFSSResponseFileHandlerCore::ReadHFSSFile()
    {
        return true;
    }
   
    double HFSSResponseFileHandlerCore::ConvolveWithFIRFilter(std::deque<double> inputBuffer)
    {
        double convolution=0.0;
        if(fFIRNBins<=0)
        {
            LERROR(lmclog,"Number of bins in the filter should be positive");
        }
        int firBinNumber=0;
        for (auto it = inputBuffer.begin();it!=inputBuffer.end(); ++it)
        {
        	convolution+=*(it)*fFilter[firBinNumber];
        	firBinNumber++;
        }

        return convolution;
    }
    
    std::pair<double,double> HFSSResponseFileHandlerCore::ConvolveWithComplexFIRFilter(std::deque<double> inputBuffer)
    {
        double convolutionMag = 0.0;
        double convolutionValueReal = 0.0;
        double convolutionValueImag = 0.0;

        if(fFIRNBins<=0)
        {
            LERROR(lmclog,"Number of bins in the filter should be positive");
        }
        int firBinNumber=0;

        for (auto it = inputBuffer.begin();it!=inputBuffer.end(); ++it)
        {
        	convolutionValueReal += *(it)*fFilterComplex[firBinNumber][0];
        	convolutionValueImag += *(it)*fFilterComplex[firBinNumber][1];
        	firBinNumber++;
        }

        std::pair<double,double> complexConvolution;
        double complexPhase = atan(convolutionValueImag/convolutionValueReal);
        double complexMag = pow(convolutionValueReal*convolutionValueReal + convolutionValueImag*convolutionValueImag, 0.5);

        /* Note that the convolution above uses a real signal to drive both the real and
         * imaginary filter components, without a phase shift.  For sinusoidal signals, this
         * seems like a reasonable approach.  It also increases the effective gain of
         * the complex filter, compared to a real filter.  To align the complex gain with
         * that of a real filter, a gain factor of 0.5x for sinusoidal drive signals is
         * applied below.  For non-sinusoidal drive signals, this scaling should not necessarily
         * be applicable.
         */

        double gainFactor = 0.5;

        return std::make_pair(gainFactor*complexMag, complexPhase);
        return complexConvolution;
    }


    TFFileHandlerCore::TFFileHandlerCore():HFSSResponseFileHandlerCore(),
    fTFComplex(NULL),
    fFIRComplex(NULL),
    fInitialTFIndex(0.0),
    fTFBinWidth(100e6)
    {
    }
    
    TFFileHandlerCore::~TFFileHandlerCore()
    {
        if (fTFComplex != NULL)
        {
            fftw_free(fTFComplex);
            fTFComplex = NULL;
        }
        
        if (fFIRComplex != NULL)
        {
            fftw_free(fFIRComplex);
            fFIRComplex = NULL;
        }
    }
    
    bool TFFileHandlerCore::Configure(const scarab::param_node& aParam)
    {
        return true;
    }


    bool TFFileHandlerCore::ConvertTFtoFIR(std::vector<std::complex<double>> &tfArray, bool GeneratedTF)
    {
        if(fTFNBins<=0)
        {
            LERROR(lmclog,"The size of transfer function has to be positive integer");
            return false;
        }
        //Might need to be moved to a different function
        fTFComplex=(fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fTFNBins);
        
        for (int i = 0; i < fTFNBins; ++i)
        {
            fTFComplex[i][0]=tfArray.at(i).real();
            fTFComplex[i][1]=tfArray.at(i).imag();
        }
        // this length is somewhat arbitrary, but works for now
        // future improvement: make it adjust depending on length of FIR

        fFIRNBins=fTFNBins + 2*fComplexFFT.GetShiftNBins();


	if(GeneratedTF)
    { 
		// If TF generated based on config file (frequency ranges from 0.9 - 1.1 times the center given in .json config file), calculate the TF bin width given number of bins (also set in .json config file)).
		double AnalyticTFBinWidth = 2./9.*fInitialTFIndex/(1.0*fTFNBins);
		fComplexFFT.SetupIFFTWindow(fTFNBins,fInitialTFIndex,AnalyticTFBinWidth, fWindowName, fWindowParam);//Uses binwidth as calculated in the previous line based on internally generated TF
	}
	else
    {
		// If TF read from externally generated TF file, use TF bin width as given in .json config file
		fComplexFFT.SetupIFFTWindow(fTFNBins,fInitialTFIndex,fTFBinWidth, fWindowName, fWindowParam); 
	}

    fFIRComplex=(fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fFIRNBins);

    fComplexFFT.ApplyWindowFunction(fTFNBins, fTFComplex);

    fComplexFFT.GenerateFIR(fTFNBins,fTFComplex,fFIRComplex);

    fResolution=fComplexFFT.GetTimeResolution();
    for (int i = 0; i < fFIRNBins; i++)
    {
        fFilter.push_back(fFIRComplex[i][0]);
    }
    fFilterComplex = fFIRComplex;

    LDEBUG( lmclog, "Finished IFFT to convert transfer function to FIR");
    return true;
    }
    
    bool TFFileHandlerCore::ReadHFSSFile()
    {
    	if(fIsFIRCreated)
    	{
    		return true;
    	}
        fTFNBins=0;
        if(!ends_with(fHFSSFilename,".txt"))
        {
            LERROR(lmclog,"The TF file " << fHFSSFilename.c_str() <<"doesn't end in .txt");
            return false;
        }
        double tfIndex;
        double tfRealValue;
        double tfImaginaryValue;
        std::vector<std::complex<double>> tfArray;
        std::fstream tfFile(fHFSSFilename.c_str(),std::ios::in);
	if (tfFile.fail()) 
	{
            LERROR(lmclog,"The TF file " << fHFSSFilename.c_str() <<" doesn't exist");
            return false;
	}
        //logic copied from /LMCPatchSignalGenerator.cc
        int totalcount=0;
        
        while (!tfFile.eof()){
            std::string lineContent;
            while(std::getline(tfFile,lineContent))
            {
                if (lineContent.find('#')==std::string::npos)
                {
                    if (totalcount%fNSkips!=0){
                        ++totalcount;
                        continue;
                    }
                    std::string token;
                    std::stringstream ss(lineContent);
                    int wordCount=0;
                    while (ss >> token)
                    {
                        if(wordCount==0)tfIndex=std::stod(token);
                        else if(wordCount==1)tfRealValue=std::stod(token);
                        else if(wordCount==2)tfImaginaryValue=std::stod(token);
                        else
                        {
                            LERROR(lmclog, "There are more column than expected in the input TF file");
                            return false;
                        }
                        ++wordCount;
                    }
                    // The TF values from HFSS are in GHz, so need to convert to Hz
                    if(fTFNBins==0)fInitialTFIndex=tfIndex*pow(10.0,9);
                    const std::complex<double> temp(tfRealValue,tfImaginaryValue);
                    tfArray.push_back(temp);
                    ++fTFNBins;
                }
            }
        }
        tfFile.close();
        LDEBUG( lmclog, "Finished reading transfer function file");
        if (!ConvertStoZ( tfArray, fConvertStoZ ) )
        {
            LERROR( lmclog, "The S parameters were not successfully converted to impedance.");
        	return false;
        }

        if (!ConvertTFtoFIR( tfArray, false ) ) //bool determines if TF was generated dynamically
        {
            return false;
        }

        if (!CropFIR( fFIRComplex, fConvertStoZ ) )
        {
            return false;
        }

        fIsFIRCreated=true;

        if (fPrintFIR)
        {

            PrintFIR( fFIRComplex, fFIRNBins, fOutputPath+"/FIRhisto.root");
            PrintFIR( fTFComplex, fTFNBins, fOutputPath+"/TFhisto.root");

            LPROG( lmclog, "Finished writing histos to output/FIRhisto.root and output/TFhisto.root");
            LPROG( lmclog, "Press Return to continue, or Cntrl-C to quit.");
            getchar();
        }

        return true;
    }

    bool TFFileHandlerCore::CropFIR(fftw_complex* anArray, bool bConvert)
    {
     	if (bConvert)
     	{
     	    int cropRange = fCropIndex;

     	    fftw_complex* aSubArray;
     	    aSubArray = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * cropRange);

            for (unsigned i=0; i<cropRange; i++)
            {
            	aSubArray[i][0] = anArray[i][0];
            	aSubArray[i][1] = anArray[i][1];
            }
            anArray = aSubArray;
            fFIRNBins = cropRange;
     	}
     	return true;
     }


    bool TFFileHandlerCore::ConvertStoZ(std::vector<std::complex<double>> &tfArray, bool bConvert)
    {
    	if (bConvert)
    	{
    	    int count = tfArray.size();
    	    double maxZ = 0.;
    	    double maxZFrequency = 0.;
    	    double qInferred = -99.;
    	    double tFrequency = fInitialTFIndex;

    	    for (int i=0; i<count; i++)
    	    {
    	        double R = tfArray.at(i).real();
    	        double X = tfArray.at(i).imag();
    	        double tfI = (1.-R*R-X*X)/((1.-R)*(1.-R)+X*X);
    	        double tfQ = -(2.*X)/((1.-R)*(1.-R)+X*X);
    	        double tfMagSquared = tfI*tfI + tfQ*tfQ;
    	        if ( (tfMagSquared > maxZ) && (tFrequency > fInitialTFIndex) )
    	        {
    	            maxZFrequency = tFrequency;
    	            maxZ = tfMagSquared;
    	            qInferred = 0.;
    	        }
    	        else if ((tfMagSquared < 0.5*maxZ) && (qInferred == 0.))
    	        {
    	            double dFrequency = tFrequency - maxZFrequency;
    	            qInferred = maxZFrequency /  (2.* dFrequency);
    	        }

    	        tfArray[i] = std::complex<double>(tfI, tfQ);
    	        tFrequency += fTFBinWidth*1.e9;
    	    }
    	    for (int i=0; i<count; i++)
    	    {
    		    tfArray[i] /= pow(maxZ,0.5);
    		    tfArray[i] *= fCharacteristicImpedance;
    	    }
    	    LPROG(lmclog, "Resonant frequency is " << maxZFrequency << " and inferred Q is " << qInferred);
    	}

    	return true;
    }

    bool TFFileHandlerCore::ConvertAnalyticGFtoFIR(std::vector<std::pair<double,std::pair<double,double> > > gfArray)
    {

    	if(fIsFIRCreated)
        {
            return true;
        }

        fFIRNBins = gfArray.size();
        fResolution = gfArray[0].first;

        fFilterComplex=(fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fFIRNBins);

        for (int i = 0; i < fFIRNBins; i++)
        {
        	fFilterComplex[i][0] = gfArray[i].second.first;
        	fFilterComplex[i][1] = gfArray[i].second.second;
        }

        if (fPrintFIR)
        {

            PrintFIR( fFilterComplex, fFIRNBins, fOutputPath+"/FIRhisto.root");
            LPROG( lmclog, "Finished writing histos to output/FIRhisto.root");
            LPROG( lmclog, "Press Return to continue, or Cntrl-C to quit.");
            getchar();
        }

        LDEBUG( lmclog, "Finished populating FIR filter with Green's function.");

    	return true;
    }


    bool TFFileHandlerCore::ConvertAnalyticTFtoFIR(double initialFreq, std::vector<std::complex<double>> tfArray)
    {

	//Replaces ReadHFSSFile() in the case where a Transfer Funciton has been generated analytically

        if(fIsFIRCreated)
        {
            return true;
        }
        fTFNBins=0;

	    fTFNBins=tfArray.size();
	    fInitialTFIndex = initialFreq;
	    if(!ConvertTFtoFIR(tfArray, true)){ //bool determines if TF was generated dynamically
            return false;
        }
        fIsFIRCreated=true;
        return true;
    }    

	bool HFSSResponseFileHandlerCore::WriteRootHisto( std::vector<double> aFilter, int nBins, bool bIQ )
	{
#ifdef ROOT_FOUND
	    char fbuffer[60];
	    if (!bIQ)
	    {
	    	int a = sprintf(fbuffer, "Real");
	    }
	    else
	    {
	    	int a = sprintf(fbuffer, "Imag");
	    }
		fRootHistoWriter->OpenFile("UPDATE");
		const char *hName = fbuffer;
		TH1D* aHisto = new TH1D(hName, "Coefficients; index; coefficient", nBins, 0., (double)nBins);
		aHisto->SetDirectory(0);

		for (unsigned i=0; i<nBins; i++)
		{
			aHisto->SetBinContent(i+1, aFilter[i]);
		}

		fRootHistoWriter->Write1DHisto(aHisto);
		fRootHistoWriter->CloseFile();
		delete aHisto;
#endif
		return true;
	}

    
    void HFSSResponseFileHandlerCore::PrintFIR( std::vector<double> aFilter, int nBins, std::string filename )
    {
    	LDEBUG( lmclog, "Printing coefficients to file ... ");
    	std::string textFile(filename);
    	std::string suffix(".txt");
    	textFile.replace(textFile.find(".root"),5,suffix);
    	FILE * fFIRout = fopen(textFile.c_str(), "w");
    	for (int i = 0; i < nBins; i++)
        {
    		fprintf(fFIRout,"%g\n", aFilter[i]);
        }
        fclose(fFIRout);
#ifdef ROOT_FOUND
        WriteRootHisto( aFilter, fFIRNBins, 0 );
#endif

    }

    void HFSSResponseFileHandlerCore::PrintFIR( fftw_complex* aFilter, int nBins, std::string filename )
    {
    	std::vector<double> vecFilter0;
    	std::vector<double> vecFilter1;
    	LDEBUG( lmclog, "Printing coefficients to file ... ");
    	std::string textFile(filename);
    	std::string suffix(".txt");
    	textFile.replace(textFile.find(".root"),5,suffix);
    	FILE * fFIRout = fopen(textFile.c_str(), "w");
    	for (int i = 0; i < nBins; i++)
    	{
    		fprintf(fFIRout,"%g %g\n", aFilter[i][0], aFilter[i][1]);
    		vecFilter0.push_back(aFilter[i][0]);
    		vecFilter1.push_back(aFilter[i][1]);
    	}
    	fclose(fFIRout);
#ifdef ROOT_FOUND
    	fRootHistoWriter = RootHistoWriter::get_instance();
    	fRootHistoWriter->SetFilename(filename);
    	fRootHistoWriter->OpenFile("RECREATE");
    	fRootHistoWriter->CloseFile();
    	WriteRootHisto( vecFilter0, nBins, 0 );
    	WriteRootHisto( vecFilter1, nBins, 1 );
#endif
    }

    FIRFileHandlerCore::FIRFileHandlerCore():HFSSResponseFileHandlerCore()
    {
    }
    
    FIRFileHandlerCore::~FIRFileHandlerCore()
    {
    }
    
    bool FIRFileHandlerCore::Configure(const scarab::param_node& aParam)
    {
        return true;
    }
    
    bool FIRFileHandlerCore::ReadHFSSFile()
    {
	if(fIsFIRCreated) 
	{
	    return true;	
	}
        fFIRNBins=0;
        if(!ends_with(fHFSSFilename,".txt"))
        {
            LERROR(lmclog,"The FIR file " << fHFSSFilename.c_str() <<" doesn't end in .txt");
            return false;
        }
        double firIndex;
        double filterMagnitude;
        FILE *firFile;
        firFile=fopen(fHFSSFilename.c_str(),"r");

	if(!access(fHFSSFilename.c_str(), F_OK ))
	{
            LERROR(lmclog,"The FIR file " << fHFSSFilename.c_str() <<" doesn't exist");
            return false;
	}
        //logic copied from /LMCPatchSignalGenerator.cc
        int count=0;
        
        while (!feof(firFile)){
            fscanf(firFile,"%lf %lf",&firIndex,&filterMagnitude);
            if (count%fNSkips==0)
            {
                fFilter.push_back(filterMagnitude);
                ++fFIRNBins;
            }
            ++count;
        }
        fclose(firFile);
        LDEBUG( lmclog, "Finished reading FIR file");
        fIsFIRCreated=true;
        return true;
    }
    
} /* namespace locust */
