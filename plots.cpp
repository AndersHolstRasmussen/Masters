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
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    for(int i = 0; i < df->GetColumnNames().size(); i++){
        cout << df->GetColumnNames()[i] << endl;
    }

    auto data = df->Define("x", "cosAngAll._cosAngAll");
    auto h = data.Histo1D({"Stats", "cos(a)", 300, -1, 1}, "x");
    auto xaxis = h->GetXaxis();
    xaxis->SetTitle("cos(\\theta)");
    xaxis->CenterTitle();
    c->SetLogy();
    h->DrawClone();
    c->Modified();
    c->Update();
    c->SaveAs("/home/anders/i257/figures/cosang.pdf");
    // c->WaitPrimitive();
    // c->Close();
}

void betaAlphaAngle(RDataFrame *df){
    auto c = new TCanvas();
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);

    int bins = 500;
    double xmin = -1;
    double xmax = 1;


    auto d0 = df->Define("d0", "betaAlphaAngle0._betaAlphaAngle0");
    auto h0 = d0.Histo1D({"Stats", "data / efficiency", bins, xmin, xmax}, "d0");

    auto d1 = df->Define("d1", "betaAlphaAngle1._betaAlphaAngle1");
    auto h1 = d1.Histo1D({"Stats", "\\beta-\\alpha angle", bins, xmin, xmax}, "d1");

    auto eff = new TH1F("stats", "Angle efficiency", bins, xmin, xmax);
    auto eff2 = new TH1F("stats", "Angle efficiency", bins, xmin, xmax);
    h0->Add(&h1.GetValue()); 

    ifstream ifile("/home/anders/i257/build/efficiencyOutput.csv");
    if (!ifile.is_open()) {
        std::cerr << "There was a problem opening the input file!\n";
        exit(1);//exit or do additional error checking
    }

    double num = 0.0;
    //keep storing values from the text file so long as data exists:
    // while (ifile >> num) {
    //     eff->Fill(num);
    // }


    while (!ifile.eof())
    {
        double angle, weight;
        ifile >> angle >> weight;
        // cout << angle << "\t " << weight << endl;
        eff->Fill(angle);
        eff2->Fill(angle, weight);
    }
    

    Double_t factor = 1.;
    eff->Scale(factor/eff->Integral(), "width");
    eff2->Scale(factor/eff2->Integral(), "width");
    h0->Scale(factor/h0->Integral(), "width"); 

    // h0->Divide(eff);
    auto xaxis = h0->GetXaxis();
    auto yaxis = h0->GetYaxis();
    xaxis->SetTitle("cos(a)");
    xaxis->CenterTitle();
    yaxis->SetTitle("count");
    yaxis->CenterTitle();
    eff->SetLineColor(kRed);
    eff->SetLineWidth(3);
    eff2->SetLineColor(kBlue);
    eff2->SetLineWidth(3);
    h0->SetLineColor(kGreen);
    h0->SetLineWidth(3);
    // h0->Divide(eff);
    cout << "Kolmogorov test: " << h0->KolmogorovTest(eff2) << endl;
    // h0->DrawClone("HIST");
    
    eff2->DrawClone("HIST L");
    // eff->DrawClone("HIST SAME L");
    // h0->Divide(eff);
    h0->DrawClone("HIST SAME L");
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

void betaAlphaDifferentEnergies(RDataFrame *df){
    auto c = new TCanvas();
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);

    int bins = 100;
    double xmin = -1;
    double xmax = 1;


    auto d0_1k = df->Define("d0", "betaAlphaAngle0._betaAlphaAngle0").Filter("E._E[0] < 1000");
    auto h0_1k = d0_1k.Histo1D({"Stats", "h0_1k", bins, xmin, xmax}, "d0");
    auto d1_1k = df->Define("d1", "betaAlphaAngle1._betaAlphaAngle1").Filter("E._E[0] < 1000");
    auto h1_1k = d1_1k.Histo1D({"Stats", "h1_1k", bins, xmin, xmax}, "d1");
    h0_1k->Add(&h1_1k.GetValue()); 


    auto d0_2k = df->Define("d0", "betaAlphaAngle0._betaAlphaAngle0").Filter("E._E[0] > 1000 && E._E[0] < 2000");
    auto h0_2k = d0_2k.Histo1D({"Stats", "h0_2k", bins, xmin, xmax}, "d0");
    auto d1_2k = df->Define("d1", "betaAlphaAngle1._betaAlphaAngle1").Filter("E._E[0] > 1000 && E._E[0] < 2000");
    auto h1_2k = d1_2k.Histo1D({"Stats", "h1_2k", bins, xmin, xmax}, "d1");
    h0_2k->Add(&h1_2k.GetValue()); 


    auto d0_3k = df->Define("d0", "betaAlphaAngle0._betaAlphaAngle0").Filter("E._E[0] > 2000 && E._E[0] < 3000");
    auto h0_3k = d0_3k.Histo1D({"Stats", "h0_3k", bins, xmin, xmax}, "d0");
    auto d1_3k = df->Define("d1", "betaAlphaAngle1._betaAlphaAngle1").Filter("E._E[0] > 2000 && E._E[0] < 3000");
    auto h1_3k = d1_3k.Histo1D({"Stats", "h1_3k", bins, xmin, xmax}, "d1");
    h0_3k->Add(&h1_3k.GetValue()); 


    auto d0_4k = df->Define("d0", "betaAlphaAngle0._betaAlphaAngle0").Filter("E._E[0] > 3000 && E._E[0] < 4000");
    auto h0_4k = d0_4k.Histo1D({"Stats", "h0_4k", bins, xmin, xmax}, "d0");
    auto d1_4k = df->Define("d1", "betaAlphaAngle1._betaAlphaAngle1").Filter("E._E[0] > 3000 && E._E[0] < 4000");
    auto h1_4k = d1_4k.Histo1D({"Stats", "h1_4k", bins, xmin, xmax}, "d1");
    h0_4k->Add(&h1_4k.GetValue()); 


    auto d0_5k = df->Define("d0", "betaAlphaAngle0._betaAlphaAngle0").Filter("E._E[0] > 4000 && E._E[0] < 5000");
    auto h0_5k = d0_5k.Histo1D({"Stats", "h0_5k", bins, xmin, xmax}, "d0");
    auto d1_5k = df->Define("d1", "betaAlphaAngle1._betaAlphaAngle1").Filter("E._E[0] > 4000 && E._E[0] < 5000");
    auto h1_5k = d1_5k.Histo1D({"Stats", "h1_5k", bins, xmin, xmax}, "d1");
    h0_5k->Add(&h1_5k.GetValue()); 


    auto d0_6k = df->Define("d0", "betaAlphaAngle0._betaAlphaAngle0").Filter("E._E[0] > 5000 && E._E[0] < 6000");
    auto h0_6k = d0_6k.Histo1D({"Stats", "h0_6k", bins, xmin, xmax}, "d0");
    auto d1_6k = df->Define("d1", "betaAlphaAngle1._betaAlphaAngle1").Filter("E._E[0] > 5000 && E._E[0] < 6000");
    auto h1_6k = d1_6k.Histo1D({"Stats", "h1_6k", bins, xmin, xmax}, "d1");
    h0_6k->Add(&h1_6k.GetValue()); 

    
    auto d0_7k = df->Define("d0", "betaAlphaAngle0._betaAlphaAngle0").Filter("E._E[0] > 6000 && E._E[0] < 7000");
    auto h0_7k = d0_7k.Histo1D({"Stats", "h0_7k", bins, xmin, xmax}, "d0");
    auto d1_7k = df->Define("d1", "betaAlphaAngle1._betaAlphaAngle1").Filter("E._E[0] > 6000 && E._E[0] < 7000");
    auto h1_7k = d1_7k.Histo1D({"Stats", "h1_7k", bins, xmin, xmax}, "d1");
    h0_7k->Add(&h1_7k.GetValue()); 

    auto d0_8k = df->Define("d0", "betaAlphaAngle0._betaAlphaAngle0").Filter("E._E[0] > 7000 && E._E[0] < 8000");
    auto h0_8k = d0_8k.Histo1D({"Stats", "h0_8k", bins, xmin, xmax}, "d0");
    auto d1_8k = df->Define("d1", "betaAlphaAngle1._betaAlphaAngle1").Filter("E._E[0] > 7000 && E._E[0] < 8000");
    auto h1_8k = d1_8k.Histo1D({"Stats", "h1_8k", bins, xmin, xmax}, "d1");
    h0_8k->Add(&h1_8k.GetValue()); 


    auto d0_2_6k = df->Define("d0", "betaAlphaAngle0._betaAlphaAngle0").Filter("E._E[0] > 1000 && E._E[0] < 6000");
    auto h0_2_6k = d0_2_6k.Histo1D({"Stats", "h0_2_6k", bins, xmin, xmax}, "d0");
    auto d1_2_6k = df->Define("d1", "betaAlphaAngle1._betaAlphaAngle1").Filter("E._E[0] > 1000 && E._E[0] < 6000");
    auto h1_2_6k = d1_2_6k.Histo1D({"Stats", "h1_2_6k", bins, xmin, xmax}, "d1");
    h0_2_6k->Add(&h1_2_6k.GetValue()); 

    auto eff = new TH1F("stats", "Angle efficiency", bins, xmin, xmax);
    ifstream ifile("/home/anders/i257/build/efficiencyOutput.csv");
    double num = 0.0;
    while (ifile >> num) {
        eff->Fill(num);
    }
    


    Double_t factor = 1.;
    eff->Scale(factor/eff->Integral(), "width");
    h0_1k->Scale(factor/h0_1k->Integral(), "width"); 
    h0_2k->Scale(factor/h0_2k->Integral(), "width");
    h0_3k->Scale(factor/h0_3k->Integral(), "width");
    h0_4k->Scale(factor/h0_4k->Integral(), "width");
    h0_5k->Scale(factor/h0_5k->Integral(), "width");
    h0_6k->Scale(factor/h0_6k->Integral(), "width");
    h0_7k->Scale(factor/h0_7k->Integral(), "width");
    h0_8k->Scale(factor/h0_8k->Integral(), "width");
    h0_2_6k->Scale(factor/h0_2_6k->Integral(), "width");
    eff->SetLineColor(kRed);
    eff->SetLineWidth(2);
    h0_1k->SetLineColor(kGreen);
    h0_2k->SetLineColor(kMagenta);
    h0_3k->SetLineColor(kBlue);
    h0_4k->SetLineColor(kYellow);
    h0_5k->SetLineColor(kBlack);
    h0_6k->SetLineColor(kPink);
    h0_7k->SetLineColor(kGreen);
    h0_8k->SetLineColor(kBlue);

    h0_2_6k->SetLineColor(kBlue);

    cout << "Kolmogorov test 0->1000: " << h0_1k->KolmogorovTest(eff) << endl;
    cout << "Kolmogorov test 1000->2000: " << h0_2k->KolmogorovTest(eff) << endl;
    cout << "Kolmogorov test 2000->3000: " << h0_3k->KolmogorovTest(eff) << endl;
    cout << "Kolmogorov test 3000->4000: " << h0_4k->KolmogorovTest(eff) << endl;
    cout << "Kolmogorov test 4000->5000: " << h0_5k->KolmogorovTest(eff) << endl;
    cout << "Kolmogorov test 5000->6000: " << h0_6k->KolmogorovTest(eff) << endl;
    cout << "Kolmogorov test 6000->7000: " << h0_7k->KolmogorovTest(eff) << endl;
    cout << "Kolmogorov test 7000->8000: " << h0_8k->KolmogorovTest(eff) << endl;
    cout << "Kolmogorov test 1000->6000: " << h0_2_6k->KolmogorovTest(eff) << endl;
    

    // h0_1k->DrawClone("HIST");
    // h0_2k->DrawClone("HIST SAME");
    // h0_3k->DrawClone("HIST SAME");
    // h0_4k->DrawClone("HIST SAME");
    // h0_5k->DrawClone("HIST SAME");
    // h0_6k->DrawClone("HIST SAME");
    h0_7k->DrawClone("Hist ");
    h0_8k->DrawClone("HIST SAME");
    // h0_2_6k->DrawClone("HIST");
    eff->DrawClone("HIST SAME");
    
}

void alphaEffeciency(RDataFrame *df){
    auto c = new TCanvas();
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);

    int bins = 100;
    double xmin = -1;
    double xmax = 1;

    auto d0 = df->Define("d0", "angToBeam._angToBeam[0]");
    auto h0 = d0.Histo1D({"Stats", "data / efficiency", bins, xmin, xmax}, "d0");

    auto d1 = df->Define("d1", "angToBeam._angToBeam[1]");
    auto h1 = d1.Histo1D({"Stats", "\\beta-\\alpha angle", bins, xmin, xmax}, "d1");

    auto eff = new TH1F("stats", "Angle efficiency", bins, xmin, xmax);
    h0->Add(&h1.GetValue()); 

    ifstream ifile("/home/anders/i257/build/efficiencyOutputAlphas.csv");
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

    eff->SetLineColor(kRed);
    h0->SetLineColor(kGreen);

    h0->DrawClone("HIST");
    eff->DrawClone("HIST SAME");

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
    TString filename = "/home/anders/i257/data/Li8/225_N102mlio.root";
    //const char *input_file = argv[2];
    chain.Add(filename);
    RDataFrame df(chain);

    // TFile *file = TFile::Open("/home/anders/i257/data/Li8/225_03N10mlio.root", "READ");
    // TTree *t; 
    // file->GetObject("tvec", t);

    // start a ROOT application window such that the plots can actually be shown
    TApplication *app = new TApplication("ROOT window", 0, 0);
    // Call the drawing of spectrum method. Default detectorId = 0

    // alphaEffeciency(&df);
    // cosang(&df);
    // EEfigure(&df);
    // betaAlphaDifferentEnergies(&df);
    betaAlphaAngle(&df);
    // individualDetectorsBetaAlphaAngle(&df);
    // betaSpec(&df);
    // angEDiff(&df);
    app->Run(); // show all canvas
    return 0;
}
