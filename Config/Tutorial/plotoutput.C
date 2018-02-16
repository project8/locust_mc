void katydidwaterfall()
{
    gStyle->SetOptStat(kFALSE);
    gStyle->SetCanvasColor(0);
    //gStyle->SetOptTitle(0);
    //gStyle->SetTitleY(0.98f);
    //gStyle->SetTitleX(0.5f);
    //gStyle->SetTitleW(0.8f);
    gStyle->SetPadLeftMargin(0.1);
    gStyle->SetPadRightMargin(0.3);
    gStyle->SetPadBottomMargin(0.14);
    gStyle->SetTitleBorderSize(0);
    gStyle->SetTitleStyle(0);
    gStyle->SetTitleSize(0.05, "t");
     

    gStyle->SetOptStat(0);
    //gStyle->SetPalette(8);

    //TCanvas *c6 = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("Waterfall"); 
    //if (c6) delete c6;
    TCanvas* c6 = new TCanvas("Waterfall", "Waterfall", 0, 0, 800, 600);
    c6->SetBorderMode(0);
    c6->SetFrameBorderMode(0);
    c6->Update();

    TFile* f = new TFile("/home/XXX/locust_mc/cbuild/bin/katydidwaterfall.root");

    int hindex = 0;
    int n;
    char buffer[30];
    //n=sprintf(buffer, "FSFFTWSpectrogram_%d", hindex);
    //n=sprintf(buffer, "PowerSpectrogram_%d", hindex);
    n=sprintf(buffer, "PSDSpectrogram_%d", hindex);
    const char *hname = buffer;
    TH1F *hwaterfall = (TH1F*)f->Get(hname);
    //hwaterfall->GetZaxis()->SetRangeUser(-0.01e-6, 0.6e-6);
    //hwaterfall->GetYaxis()->SetRangeUser(84.e6,86.e6);
    hwaterfall->GetXaxis()->SetRangeUser(0.0,0.001);
    hwaterfall->SetNdivisions(8, "x");
    //hwaterfall->GetYaxis()->SetRangeUser(64.5e6,65.e6);
    //hwaterfall->GetZaxis()->SetRangeUser(-0.1e-23, 6.5e-21);
    hwaterfall->GetZaxis()->SetRangeUser(0., 1.e-22);

    hwaterfall->Scale(50.);  // take out 50 ohm resistor built into katydid.


    hwaterfall->GetYaxis()->SetTitleSize(0.05);
    hwaterfall->GetXaxis()->SetTitleSize(0.05);
    hwaterfall->SetTitle("spectrogram");
    hwaterfall->GetYaxis()->SetTitle("frequency");
    hwaterfall->GetZaxis()->SetTitleOffset(1.2);
    //hwaterfall->GetZaxis()->SetTitle("volts");
    //hwaterfall->SetTitle("30 keV electron in harmonic trap, #theta=90.0^{o}, mixed down to -50 MHz");
    //hwaterfall->SetTitle("89^{o} < #theta < 90^{o}, noise power 1.e-25 W/Hz");


    //threshold cut
    // 
    // for (int i=0; i<hwaterfall->GetXaxis()->GetNbins(); i++)
    //   for (int j=0; j<hwaterfall->GetYaxis()->GetNbins(); j++)
    //     if (hwaterfall->GetBinContent(i+1, j+1) < 1.e-22)
    //       hwaterfall->SetBinContent(i+1, j+1, 0.);
    // 


    hwaterfall->DrawCopy("colz");

    c6->Update();
    hwaterfall->DrawCopy("colz");
   c6->Update();

}



TGraph* BuildKassGraphGuidingCenterXYPlane(const char* filename)
{

    TFile* f = new TFile(string(filename).c_str());

    double* XArray = new double[10000000];
    double* YArray = new double[10000000];
    int nentries = 0;
    double time = 0.;
    int counter = 0;


    t1 = (TTree *)f->FindObjectAny("component_step_world_DATA");
    if (!t1) {printf("can't find the tree\n");}
    else
    {
        nentries = (Int_t)t1->GetEntries();
        for (Int_t i=0; i<nentries; i++)
        {
            t1->GetEntry(i);
            if ((i%1000==0))
            {
                YArray[counter] = t1->GetLeaf("guiding_center_position_y")->GetValue(0);
                XArray[counter] = t1->GetLeaf("guiding_center_position_x")->GetValue(0);
                counter += 1;
                time = t1->GetLeaf("time")->GetValue(0);
            }
        }
    } // else                                                                   

    printf("counter is %d and old_time is %f\n", counter, time);
    TGraph *grZ = new TGraph(counter, XArray, YArray);

    return grZ;

}


TGraph* BuildKassGraphXYPlane(const char* filename)
{

    TFile* f = new TFile(string(filename).c_str());

    double* XArray = new double[10000000];
    double* YArray = new double[10000000];
    int nentries = 0;
    double time = 0.;
    int counter = 0;

    t1 = (TTree *)f->FindObjectAny("component_step_world_DATA");
    if (!t1) {printf("can't find the tree\n");}
    else  
    {
        nentries = (Int_t)t1->GetEntries();
        for (Int_t i=0; i<nentries; i++) 
        {
            t1->GetEntry(i);
            if ((i%10001==0))
            {        
                YArray[counter] = t1->GetLeaf("position_y")->GetValue(0);
                XArray[counter] = t1->GetLeaf("position_x")->GetValue(0);
                counter += 1;
                time = t1->GetLeaf("time")->GetValue(0);

            }
        }
    } // else

    printf("counter is %d and old_time is %f\n", counter, time);
    TGraph *grZ = new TGraph(counter, XArray, YArray);

    return grZ;
}


void KassPlotXYPlane()
{

    gStyle->SetOptStat(kFALSE);
    gStyle->SetCanvasColor(0);
    //gStyle->SetOptTitle(0);
    gStyle->SetTitleY(0.98f);
    gStyle->SetTitleX(0.5f);
    gStyle->SetTitleW(0.8f);
    gStyle->SetTitleBorderSize(0);
    gStyle->SetTitleStyle(0);
    //gStyle->SetTitleSize(0.06, "t");
    gStyle->SetTitleX(0.55f);
    gStyle->SetPadLeftMargin(0.22); 

    TCanvas* c133 = new TCanvas("Kassiopeia phi", "Kassiopeia phi", 0, 0, 800, 600);
    c133->SetBorderMode(0);
    c133->SetFrameBorderMode(0);
    c133->Update();

    TLine* line1 = new TLine(0., 0., 0., 0.);
    TH2D *hXY = new TH2D("hXY", "Trajectory in X Y plane", 500, -0.006, 0.006, 500, -0.006, 0.006);
    hXY->GetYaxis()->SetTitleOffset(1.5);
    // hXY->GetXaxis()->SetRangeUser(-0.00125,-0.00075);
    // hXY->GetYaxis()->SetRangeUser(-0.0001,0.0001);
    hXY->DrawCopy();


    string filename[1];
    filename[0] = "/home/XXX/locust_mc/cbuild/output/Kassiopeia/Phase2Output.root";


    TLatex* t = new TLatex(0., 0., "test");

    for (int iter=0; iter<1; iter++)
    {

        printf("iter is %d\n", iter);

        const char *file = filename[iter].c_str();

        TGraph *grZ = BuildKassGraphXYPlane(file);
        TGraph *grGuidingCenter = BuildKassGraphGuidingCenterXYPlane(file);
        grZ->GetXaxis()->SetTitle("X (m)");
        grZ->GetYaxis()->SetTitleOffset(1.3);
        grZ->GetXaxis()->SetTitleOffset(1.0);
        grZ->GetYaxis()->SetTitleSize(0.05);
        grZ->GetXaxis()->SetTitleSize(0.05);
        grZ->GetYaxis()->SetLabelSize(0.045);
        grZ->GetXaxis()->SetLabelSize(0.045);
        grZ->GetYaxis()->SetTitle("Y (m)");
        grZ->SetTitle("XY plane");
        //grZ->GetXaxis()->SetRangeUser(-0.006, 0.006);
        //grZ->GetYaxis()->SetRangeUser(-0.003, 0.003);
        grZ->SetLineWidth(1);
        grZ->SetLineColor(iter+1);
        grZ->SetMarkerColor(iter+1);
        if (iter==1) grZ->SetLineColor(2);
        grZ->Draw("lsame");
        grGuidingCenter->SetLineColor(3);
        grGuidingCenter->Draw("lsame");

        line1->SetLineColor(iter+1);
        line1->SetLineWidth(1);
        //line1->DrawLine(-0.002, 0.0055-(double)iter*0.0005, -0.001, 0.0055-(double)iter*0.0005);
        t->SetTextFont(42);
        t->SetTextSize(0.04);
        //t->DrawLatex(-0.004, -0.002, "#theta=89.8^{o}, 0 < t < 1.e-4 s");

    } //iter

    e = new TEllipse(0.,0.,0.00502920,0.00502920);
    e->SetFillStyle(0);
    e->SetLineWidth(3.);
    e->SetLineStyle(2);
    e->Draw();
    line1->SetLineWidth(1.);
    line1->SetLineColor(1);
    line1->SetLineStyle(2);
    line1->DrawLine(-0.0050292, 0., 0.0050292, 0.);
    line1->DrawLine(0., -0.0050292, 0., 0.0050292);
    t->SetTextFont(42);
    t->SetTextSize(0.04);
    t->DrawLatex(0.0015, -0.0051, "waveguide");


    if (0 == 1) // phase 1
    {
        line1->SetLineWidth(3.);
        line1->SetLineColor(1);
        line1->SetLineStyle(2);
        line1->DrawLine(-0.005334, 0.002159, 0.005334, 0.002159);
        line1->DrawLine(-0.005334, -0.002159, -0.005334, 0.002159);
        line1->DrawLine(0.005334, -0.002159, 0.005334, 0.002159);
        line1->DrawLine(-0.005334, -0.002159, 0.005334, -0.002159);

        t->SetTextFont(42);
        t->SetTextSize(0.03);
        t->DrawLatex(-0.005, 0.0023, "waveguide removed for large orbits");
    }

    //TH1D* hfft = FFT(ZArray, 100000);
    //hfft->DrawCopy();

     
}




    // plot Bz.  Use radial cut to pick center axis.
void fieldmap() {

    double pz_array[200000];
    double Bz_array[200000];
    double x=0., y=0.;

    TFile *f = TFile::Open("~/locust_mc/cbuild/output/Kassiopeia/Phase1Output.root");

    // Create a TTreeReader for the tree, for instance by passing the
    // TTree's name and the TDirectory / TFile it is in.
    TTreeReader myReader("component_step_world_DATA", f);
    TTreeReaderValue<double> pz(myReader,"position_z");
    TTreeReaderValue<double> Bz(myReader,"magnetic_field_z");
    TTreeReaderValue<double> px(myReader,"position_x");
    TTreeReaderValue<double> py(myReader,"position_y");

    int counter = 0;
    while (myReader.Next()) 
    {
        x = *px;
        y = *py;
        if (x*x+y*y < 0.00001)  // pick central field line
        {
            pz_array[counter] = *pz;
            Bz_array[counter] = *Bz;
            counter += 1;
        }
    }

    TGraph *grField = new TGraph(counter, pz_array, Bz_array);

    grField->SetTitle("B_{z} on axis");
    grField->GetYaxis()->SetTitle("B_{z} (T)");
    grField->GetXaxis()->SetTitle("z position (m)");
    grField->GetYaxis()->SetTitleSize(0.05);
    grField->GetXaxis()->SetTitleSize(0.05);
    grField->GetYaxis()->SetLabelSize(0.045);
    grField->GetXaxis()->SetLabelSize(0.045);

    grField->GetYaxis()->SetTitleOffset(1.2);
    grField->SetMarkerStyle(8);
    grField->Draw("ap");

    //TF1 *fa1 = new TF1("fa1","58.*x*x+0.9957",-0.006,0.006); 
    //fa1->Draw("same"); // parabola line.
}

