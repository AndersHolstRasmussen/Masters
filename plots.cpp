#include <TChain.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TFile.h>
#include <TStyle.h>
#include <TF1.h>
#include <TGraph.h>
#include <TVectorD.h>
#include <TLegend.h>
#include <TSpectrum.h>
#include <TApplication.h>
#include <string>
#include <sstream>
#include <TColor.h>
#include <TH2D.h>
#include <TRandom.h>
#include <TPaveText.h>
#include <ROOT/RDataFrame.hxx>
#include <TTree.h>
#include <iomanip>
#include <cstdio>
#include <ctime>
#include <TColorGradient.h>
#include <ROOT/RCsvDS.hxx>
#include <fstream>

using namespace std;
using namespace ROOT;


void EEfigure(RDataFrame *df){
    // auto c = new TCanvas("c", "k", 200, 110, 800, 800);
    auto c = new TCanvas();
    c->SetTitle("EE");
    c->SetGridx();
    c->SetGridy();
    gStyle->SetPalette(kViridis);
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);

    auto a1 = df->Define("a1", "E._E[0]").Define("a2", "E._E[1]");

    auto h = a1.Histo2D({"stats", "Energy vs energy", 300, 0, 8000, 300, 0, 8000}, "a1", "a2");
    auto xaxis = h->GetXaxis();
    auto yaxis = h->GetYaxis();
    auto zaxis = h->GetZaxis();
    // h->SetContour(1000);
    c->SetLogz();
    xaxis->SetTitle(" E_{\\alpha1} [keV]");
    xaxis->CenterTitle();
    yaxis->CenterTitle();
    yaxis->SetTitle(" E_{\\alpha2} [keV]");
    h->DrawClone("contz");


    c->Modified();
    c->Update();
    c->SaveAs("/home/anders/i257/figures/EE.pdf");
    // // c->WaitPrimitive();
    // // c->Close();
    
}

void cosang(RDataFrame *df){
    auto c = new TCanvas();
    for(int i = 0; i < df->GetColumnNames().size(); i++){
        cout << df->GetColumnNames()[i] << endl;
    }

    auto data = df->Define("x", "cosA._cosA");
    auto h = data.Histo1D({"Stats", "cos(a)", 300, -1.05, 1.05}, "x");
    h->DrawClone();
    c->Modified();
    c->Update();
    c->SaveAs("/home/anders/i257/figures/cosang.pdf");
    // c->WaitPrimitive();
    // c->Close();
}

void betaAlphaAngle(RDataFrame *df){
    auto c = new TCanvas();
    // gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);

    int bins = 100;
    double xmin = -1;
    double xmax = 1;


    auto d0 = df->Define("d0", "betaAlphaAngle0._betaAlphaAngle0");
    auto h0 = d0.Histo1D({"Stats", "data / efficiency", bins, xmin, xmax}, "d0");

    auto d1 = df->Define("d1", "betaAlphaAngle1._betaAlphaAngle1");
    auto h1 = d1.Histo1D({"Stats", "\\beta-\\alpha angle", bins, xmin, xmax}, "d1");

    auto eff = new TH1F("stats", "Angle efficiency", bins, xmin, xmax);
    h0->Add(&h1.GetValue()); 

    ifstream ifile("/home/anders/i257/build/efficiencyOutput.csv");
    if (!ifile.is_open()) {
        std::cerr << "There was a problem opening the input file!\n";
        exit(1);//exit or do additional error checking
    }

    double num = 0.0;
    //keep storing values from the text file so long as data exists:
    while (ifile >> num) {
        eff->Fill(num);
    }

    Double_t factor = 1.;
    eff->Scale(factor/eff->Integral(), "width");
    h0->Scale(factor/h0->Integral(), "width"); 

    // h0->Divide(eff);
    auto xaxis = h0->GetXaxis();
    auto yaxis = h0->GetYaxis();
    xaxis->SetTitle("cos(a)");
    xaxis->CenterTitle();
    yaxis->SetTitle("count");
    yaxis->CenterTitle();
    eff->SetLineColor(kRed);
    h0->SetLineColor(kGreen);

    cout << "Kolmogorov test: " << h0->KolmogorovTest(eff) << endl;
    eff->DrawClone("HIST");
    h0->DrawClone("HIST SAME");
    // h0->Divide(eff);
    // h0->DrawClone("HIST");
    // eff->DrawClone();
    c->Modified();
    c->Update();
    // c->SaveAs("/home/anders/i257/figures/dataDivEff.pdf");

}

void individualDetectorsBetaAlphaAngle(RDataFrame *df){
    auto c = new TCanvas();
    int bins = 100;
    double xmin = -1;
    double xmax = 1;
    
    int detectorNr = 1;
    string x = to_string(detectorNr);
    string title = "Angle where beta hit in detector 1 or " + to_string(detectorNr);
    auto d0 = df->Define("d0", "betaAlphaAngle0._betaAlphaAngle0").Filter("ibeta._ibeta[0] == " + x); // + "|| ibeta._ibeta[0] == 1"
    auto h0 = d0.Histo1D({"Stats", title.c_str() , bins, xmin, xmax}, "d0");

    auto d1 = df->Define("d1", "betaAlphaAngle1._betaAlphaAngle1").Filter("ibeta._ibeta[0] == " + x); //  + "|| ibeta._ibeta[0] == 1"
    auto h1 = d1.Histo1D({"Stats", "\\beta-\\alpha angle", bins, xmin, xmax}, "d1");

    auto eff = new TH1F("stats", "Angle efficiency", bins, xmin, xmax);
    ifstream ifile("/home/anders/i257/build/efficiencyOutput.csv");


    if (!ifile.is_open()) {
        std::cerr << "There was a problem opening the input file!\n";
        exit(1);//exit or do additional error checking
    }

    double num = 0.0;
    //keep storing values from the text file so long as data exists:
    while (ifile >> num) {
        eff->Fill(num);
    }
    h0->Add(&h1.GetValue());
    cout << "Kolmogorov test: " << h0->KolmogorovTest(eff) << endl;
    h0->Divide(eff);
    
    h0->DrawClone();
    // eff->DrawClone();
    // string savetitle = "/home/anders/i257/figures/betaAngles/angleWhereBetaWasInDet" + to_string(detectorNr) + ".pdf";
    // c->SaveAs(savetitle.c_str());
}

void betaSpec(RDataFrame *df){
    auto c = new TCanvas();
    auto data = df->Define("x", "EPbeta._EPbeta");
    auto h = data.Histo1D({"Stats", "Beta spectrum (from pads)", 300, 1, 2000}, "x");
    TAxis *xaxis = h->GetXaxis();
    xaxis->SetTitle("name");
    h->DrawClone();
    c->Modified();
    c->Update();
    c->SaveAs("/home/anders/i257/figures/betaSpec.pdf");
}

void angEDiff(RDataFrame *df){
    
    auto c = new TCanvas("c", "k", 200, 110, 1800, 1000);
    c->SetTitle("Energy difference and angle");

    auto cosang = df->Define("c", "cosA._cosA").Take<double>("c").GetValue();
    auto e0 = df->Define("e0", "E._E[0]").Take<double>("e0").GetValue();
    auto e1 = df->Define("e1", "E._E[1]").Take<double>("e1").GetValue();
    int n = e0.size();
    double *d;
    d = new double[n];

    for(int i = 0; i < n; i++){
        d[i] = abs(e0[i] - e1[i]);
    }
    cout << n << endl;

    TGraph *g = new TGraph (n, &cosang[0], &d[0]);
    g->SetMarkerStyle(1);
    g->SetMarkerColorAlpha(kBlack, 1);
    g->GetHistogram()->SetMaximum(8000);
    g->GetHistogram()->SetMinimum(0);
    
    g->Draw("AP");
    c->SaveAs(("/home/anders/i257/figures/angEDiff.png"));
    c->Modified();
    c->Update();

}

int main(int argc, char *argv[]) {
    ROOT::EnableImplicitMT(8);
    int detectorId;
    std::setw(2);
    std::setprecision(5);

    if (argc > 1) { // If an input has been given
        istringstream ss(argv[1]);  //Convert the input to int
        ss >> detectorId;           //Convert the input to int
    } 

    TChain chain("tree");
    TString filename = "/home/anders/i257/data/Li8/225_03N10mlio.root";
    //const char *input_file = argv[2];
    chain.Add(filename);
    RDataFrame df(chain);

    // TFile *file = TFile::Open("/home/anders/i257/data/Li8/225_03N10mlio.root", "READ");
    // TTree *t; 
    // file->GetObject("tvec", t);

    // start a ROOT application window such that the plots can actually be shown
    TApplication *app = new TApplication("ROOT window", 0, 0);
    // Call the drawing of spectrum method. Default detectorId = 0

    // cosang(&df);
    // EEfigure(&df);
    betaAlphaAngle(&df);
    // individualDetectorsBetaAlphaAngle(&df);
    // betaSpec(&df);
    // angEDiff(&df);
    app->Run(); // show all canvas
    return 0;
}
