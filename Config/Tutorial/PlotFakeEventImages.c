/*
Author: Penny Slocum
Date: 10/11/18

ROOT script to make waterfall plots from generated fake tracks:
root -l
.L PlotFakeTrackImages.c
PlotImages()
*/

TH2D* GetSpectrogram(const char* filename)
{
TFile* f = TFile::Open(filename);
TH2D *hwaterfall = (TH2D*)f->Get("PSDSpectrogram_0_0");
return hwaterfall;
}

TH2D* GetLabels(TH2D* hspectrogram, double threshold)
{
  TH2D *hlabels = (TH2D*)hspectrogram->Clone("hlabels");
  for (unsigned i=0; i<hspectrogram->GetXaxis()->GetNbins(); i++)
    {
    for (unsigned j=0; j<hspectrogram->GetYaxis()->GetNbins(); j++)
      {
	if (hspectrogram->GetBinContent(i+1,j+1)>threshold)
	            hlabels->SetBinContent(i+1,j+1, 1.0);
        else
          hlabels->SetBinContent(i+1,j+1, 0.0);            
      }
    }
  return hlabels;
}


TGraph* GetLabelGraph(TH2D* hlabels)
{
  int count=0;
  int nbinstime = hlabels->GetXaxis()->GetNbins();
  int nbinsfreq = hlabels->GetYaxis()->GetNbins();
  double* time = new double[nbinstime];
  double* freq = new double[nbinsfreq];

  for (unsigned i=0; i<hlabels->GetXaxis()->GetNbins(); i++)
    {
      for (unsigned j=0; j<hlabels->GetYaxis()->GetNbins(); j++)
	{
	  if (hlabels->GetBinContent(i+1,j+1)>0.)
            {
              time[count] = hlabels->GetXaxis()->GetBinCenter(i+1);
              freq[count] = hlabels->GetYaxis()->GetBinCenter(j+1);
	      count += 1;
	      // printf("found one %g %g\n", time[count-1], freq[count-1]);             
	    }	    	    
	}
    }
  TGraph* grlabels = new TGraph(count, time, freq);
  delete[] time;
  delete[] freq;

  return grlabels;
}


void PrintLabels(TGraph* grlabels)
{
  grlabels->Print();

}


void PlotImages()
{
gStyle->SetOptStat(0);
gStyle->SetPadRightMargin(0.2);
string filename[1];
char buffer[100];
char buffertext[100];

for (int seed=55; seed<56; seed++)
  {
    int n=sprintf(buffer, "/home/les67/locust_faketrack_waterfall.root");

    const char *file = buffer;

    TFile* f = TFile::Open(file);
    if (!((!f)||f->IsZombie()))
    {

    TH2D* hspectrogram = GetSpectrogram(file);
    TH2D* hlabels = GetLabels(hspectrogram, 40e-21);  // threshold for labeling goes here.
//    TGraph* grlabels = GetLabelGraph(hlabels);
//    PrintLabels(grlabels);  // print labels to terminal.

    hspectrogram->GetXaxis()->SetRangeUser(0.,0.02);
    hspectrogram->GetYaxis()->SetRangeUser(149.e6, 170.e6);

    TCanvas *c = new TCanvas;
    n=sprintf(buffertext, "hspectrogram_%d.png", seed);
    const char *pngname = buffertext;
    hspectrogram->GetYaxis()->SetTitleOffset(1.25);
    hspectrogram->DrawCopy("colz");

/*
    grlabels->SetMarkerColor(2);
    grlabels->SetMarkerStyle(8);
    grlabels->SetMarkerSize(0.6);
    grlabels->SetLineWidth(3.);
    grlabels->SetLineColor(2);
    grlabels->Draw("psame");
*/


 /*  
    TImage *img = TImage::Create();
    img->FromPad(c);
    img->WriteImage(pngname);
    delete img;
    delete hspectrogram;
    delete c;
*/
    f->Close();
    }

  }

}
