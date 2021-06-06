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
#include <TRatioPlot.h>

using namespace std;
using namespace ROOT;

void singleDetectorEff(RDataFrame *df){
    auto c = new TCanvas();
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kViridis);
    int bins = 100;
    double xmin = 0;
    double xmax = 1;

    auto eff = new TH1F("stats", "Angle efficiency", bins, xmin, xmax);
    auto eff2 = new TH1F("stats", "Angle efficiency", bins, xmin, xmax);

    ifstream ifile("/home/anders/i257/build/angEffDet2.csv");
    if (!ifile.is_open()) {
        std::cerr << "There was a problem opening the input file!\n";
        exit(1);//exit or do additional error checking
    }
    while (!ifile.eof())
    {
        double angle, weight;
        ifile >> angle >> weight;
        eff->Fill(angle);
        eff2->Fill(angle, weight);
    }
    // Double_t factor = 1.;
    // eff->Scale(factor/eff->Integral(), "width");
    // eff2->Scale(factor/eff2->Integral(), "width"); 

    eff->SetLineColor(kRed);
    eff->SetLineWidth(3);
    eff2->SetLineColor(kBlue);
    eff2->SetLineWidth(3);
    auto xaxis = eff2->GetXaxis();
    auto yaxis = eff2->GetYaxis();
    auto zaxis = eff2->GetZaxis();
    // h->SetContour(1000);
    c->SetLogz();
    xaxis->SetTitle("cos(\\theta)");
    xaxis->CenterTitle();
    yaxis->CenterTitle();
    yaxis->SetTitle("Count");

    eff2->DrawClone("HIST");
    // eff->DrawClone("HIST");

    c->SaveAs("/home/anders/i257/figures/det2Eff.pdf");

}

void mexiHatTheory(RDataFrame *df){
    auto c = new TCanvas("c", "k", 200, 110, 700, 700);
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kViridis);
    int bins = 16;
    double xmin = 0;
    double xmax = 16;

    auto eff = new TH1F("stats", "Angle efficiency", bins, xmin, xmax);
    auto some = new TH2F("stats", "title", bins, xmin, xmax, bins, xmin, xmax);
    ifstream ifile("/home/anders/i257/build/effDet2Corrected.csv");
    if (!ifile.is_open()) {
        std::cerr << "There was a problem opening the input file!\n";
        exit(1);//exit or do additional error checking
    }
    double w[256];
    int k = 0;
    double num = 0.0;
   
    while (!ifile.eof())
    {
        double weight;
        ifile >> weight;
        w[k] = weight;
        k++;
    }

    k=0;
    double fi[16], bi[16];
    for (int i = 0; i < 16; i++){
        for (int j = 0; j < 16; j++){
            fi[i] = i;
            bi[j] = j;
            some->Fill(fi[i], bi[j], w[k]);
            k++;
        }
    }
    

    some->DrawClone("col");
    c->Modified();
    c->Update();
    c->WaitPrimitive();
    c->Close();
    c->SaveAs("/home/anders/i257/figures/mexihatDet2Corrected.pdf");
}

void mexiHatDetector(RDataFrame *df){
    auto c = new TCanvas("c", "k", 200, 110, 700, 700);
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kViridis);
    int bins = 16;
    double xymin = 0;
    double xymax = 17;

    auto data = df->Define("x", "FI._FI[0]").Define("y", "BI._BI[0]").Filter("ibeta._ibeta[0] == 5 "); // || i._i[1] == 1
    auto h = data.Histo2D({"stats", "title", bins, xymin, xymax, bins, xymin, xymax}, "x", "y");
    auto xaxis = h->GetXaxis();
    auto yaxis = h->GetYaxis();
    auto zaxis = h->GetZaxis();
    xaxis->SetTitle("Front strip index ");
    xaxis->CenterTitle();
    yaxis->CenterTitle();
    yaxis->SetTitle(" Back strip index");
    h->DrawClone("col");
    c->Modified();
    c->Update();
    // c->SaveAs("/home/anders/i257/figures/mexihatDet2.pdf");

}

void detectorEff(RDataFrame *df){
    auto c = new TCanvas();
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    int bins = 300;
    double xmin = -1;
    double xmax = 1;

    auto eff = new TH1F("stats", "Angle efficiency", bins, xmin, xmax);
    auto eff2 = new TH1F("stats", "Angle efficiency", bins, xmin, xmax);

    ifstream ifile("/home/anders/i257/build/efficiencyOutputAllDet.csv");
    if (!ifile.is_open()) {
        std::cerr << "There was a problem opening the input file!\n";
        exit(1);//exit or do additional error checking
    }

    double num = 0.0;



    while (!ifile.eof())
    {
        double angle, weight;
        ifile >> angle >> weight;
        // cout << angle << "\t " << weight << endl;
        eff->Fill(angle);
        eff2->Fill(angle, weight);
    }
    auto xaxis = eff2->GetXaxis();
    auto yaxis = eff2->GetYaxis();
    xaxis->SetTitle("cos(a)");
    xaxis->CenterTitle();
    yaxis->SetTitle("count");
    yaxis->CenterTitle();
    eff->SetLineColor(kRed);
    eff->SetLineWidth(3);
    eff2->SetLineColor(kBlue);
    eff2->SetLineWidth(3);

    Double_t factor = 1.;
    eff->Scale(factor/eff->Integral(), "width");
    eff2->Scale(factor/eff2->Integral(), "width");

    // eff2->DrawClone("HIST");
    eff2->DrawClone("HIST");
    eff->DrawClone("HIST SAME");

   auto legend = new TLegend(0.15, 0.8, 0.5, 0.9);
   legend->AddEntry(eff,"Angular efficiency without effective pixel area", "L");
   legend->AddEntry(eff2,"Angular efficiency");
   legend->SetTextSize(0.02);
   legend->Draw();


    c->SaveAs("/home/anders/i257/figures/allDetEff.pdf");
    c->Modified();
    c->Update();


}

void EEfigure(RDataFrame *df){
    // auto c = new TCanvas("c", "k", 200, 110, 800, 800);
    auto c = new TCanvas();
    c->SetTitle("EE");
    c->SetGridx();
    c->SetGridy();
    gStyle->SetPalette(kViridis);
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    double min = 0, max = 8000;
    int bins = 300;

    auto a1 = df->Define("a1", "E._E[0]").Define("a2", "E._E[1]");

    auto h = a1.Histo2D({"stats", "Energy vs energy", bins, min, max, bins, min, max}, "a1", "a2");
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
    c->SaveAs("/home/anders/i257/figures/EEAngAndMoment.pdf");
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

    int bins = 300;
    double xmin = -1;
    double xmax = 1;


    auto d0 = df->Define("d0", "betaAlphaAngle0._betaAlphaAngle0");
    auto h0 = d0.Histo1D({"Stats", "data / efficiency", bins, xmin, xmax}, "d0");

    auto d1 = df->Define("d1", "betaAlphaAngle1._betaAlphaAngle1");
    auto h1 = d1.Histo1D({"Stats", "\\beta-\\alpha angle", bins, xmin, xmax}, "d1");

    auto eff = new TH1F("stats", "Angle efficiency", bins, xmin, xmax);
    auto eff2 = new TH1F("stats", "Angle efficiency", bins, xmin, xmax);
    h0->Add(&h1.GetValue()); 

    ifstream ifile("/home/anders/i257/build/effBetaDets.csv");
    if (!ifile.is_open()) {
        std::cerr << "There was a problem opening the input file!\n";
        exit(1);//exit or do additional error checking
    }
    while (!ifile.eof())
    {
        double angle, weight;
        ifile >> angle >> weight;
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

    // h0->Divide(eff2);



    
    cout << "Kolmogorov test: " << h0->KolmogorovTest(eff2, "N") << endl;
    eff2->DrawClone("HIST");
    h0->DrawClone("HIST SAME");
    
    auto hclone = h0->Clone();
    auto effclone = eff2->Clone();
    
    auto legend = new TLegend(0.15, 0.8, 0.3, 0.9);
    legend->AddEntry(hclone, "Measured angle");
    legend->AddEntry(eff2,"Angular efficiency");
    legend->SetTextSize(0.02);
    // legend->Draw();
    c->Modified();
    c->Update();
    // c->SaveAs("/home/anders/i257/figures/betaAngles/dataDivEffCenterCorrected.pdf");

}

void individualDetectorsBetaAlphaAngle(RDataFrame *df){
    auto c = new TCanvas();
    int bins = 100;
    double xmin = -1;
    double xmax = 1;
    
    int detectorNr = 1;
    string x = to_string(detectorNr);
    string title = "Angle where beta hit in detector 1 or " + to_string(detectorNr);
    auto d0 = df->Define("d0", "betaAlphaAngle0._betaAlphaAngle0").Filter("ibeta._ibeta[0] == 3" + x); // + "|| ibeta._ibeta[0] == 1"
    auto h0 = d0.Histo1D({"Stats", title.c_str() , bins, xmin, xmax}, "d0");

    auto d1 = df->Define("d1", "betaAlphaAngle1._betaAlphaAngle1").Filter("ibeta._ibeta[0] == 3" + x); //  + "|| ibeta._ibeta[0] == 1"
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
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kViridis);
    auto dataD = df->Define("x", "Ebeta._Ebeta").Filter("ibeta._ibeta[0] == 5 "); // || ibeta._ibeta[1] == 5 || ibeta._ibeta[2] == 5 || ibeta._ibeta[3] == 5 || ibeta._ibeta[4] == 5 || ibeta._ibeta[5] == 5 || ibeta._ibeta[6] == 5
    auto data2 = df->Define("x", "Ebeta._Ebeta").Filter("ibeta._ibeta[0] == 1 "); // || ibeta._ibeta[1] == 1 || ibeta._ibeta[2] == 1 || ibeta._ibeta[3] == 1 || ibeta._ibeta[4] == 1 || ibeta._ibeta[5] == 1 || ibeta._ibeta[6] == 1
    auto datapadD = df->Define("x", "EPbeta._EPbeta").Filter("ibeta._ibeta[0] == 5").Filter("EPbeta._EPbeta[0] > 1");
    auto datapad2 = df->Define("x", "EPbeta._EPbeta").Filter("ibeta._ibeta[0] == 1").Filter("EPbeta._EPbeta[0] > 1");
    auto hD = dataD.Histo1D({"Stats", "Beta spectrum", 300, 1, 2000}, "x");
    auto h2 = data2.Histo1D({"Stats", "Beta spectrum", 300, 1, 2000}, "x");
    auto hpadD = datapadD.Histo1D({"Stats", "Beta spectrum", 300, 1, 2000}, "x");
    auto hpad2 = datapad2.Histo1D({"Stats", "Beta spectrum", 300, 1, 2000}, "x");
    TAxis *xaxis = hD->GetXaxis();
    hD->SetLineColor(kRed);
    hD->SetLineWidth(3);
    h2->SetLineColor(kBlue);
    h2->SetLineWidth(3);
    hpadD->SetLineColor(kPink);
    hpadD->SetLineWidth(3);
    hpad2->SetLineColor(kGreen);
    hpad2->SetLineWidth(3);
    Double_t factor = 1.;
    hD->Scale(factor/hD->Integral(), "width");
    h2->Scale(factor/h2->Integral(), "width");
    // hpadD->Scale(factor/hpadD->Integral(), "width");
    hpad2->Scale(factor/hpad2->Integral(), "width");

    xaxis->SetTitle("E_{\\beta} [keV]");
    

    hpad2->DrawClone("HIST");
    hD->DrawClone("HIST SAME");
    h2->DrawClone("HIST SAME");
    // hpadD->DrawClone("HIST SAME");
    


    auto legend = new TLegend(0.78, 0.75, 0.88, 0.88);
    legend->AddEntry(hD->Clone(),"DetD");
    legend->AddEntry(h2->Clone(),"Det2");
    // legend->AddEntry(hpadD->Clone(), "PadD");
    legend->AddEntry(hpad2->Clone(), "Pad2");
    legend->SetTextSize(0.03);
    legend->Draw();

    c->Modified();
    c->Update();
    c->SaveAs("/home/anders/i257/figures/betaSpec.pdf");
}

void alphaAndBetaEnergy(RDataFrame *df){
    auto c = new TCanvas();
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kViridis);
    auto data2 = df->Define("x", "Ebeta._Ebeta").Filter("ibeta._ibeta[0] == 1 ");
    auto dataa = df->Define("x", "E._E[0]").Filter("i._i[0] == 1");
    auto dataa2 = df->Define("x", "E._E[1]").Filter("i._i[1] == 1");
    auto h2 = data2.Histo1D({"Stats", "Beta spectrum", 300, 1, 6000}, "x");
    auto ha = dataa.Histo1D({"Stats", "Beta spectrum", 300, 1, 6000}, "x");
    auto ha2 = dataa2.Histo1D({"Stats", "Beta spectrum", 300, 1, 6000}, "x");
    ha->Add(&ha2.GetValue());
    h2->SetLineColor(kBlue);
    h2->SetLineWidth(3);
    ha->SetLineColor(kRed);
    ha->SetLineWidth(3);
    Double_t factor = 1.;
    // ha->Scale(factor/ha->Integral(), "width");
    // h2->Scale(factor/h2->Integral(), "width");
    h2->DrawClone("hist");
    ha->DrawClone("hist same");
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

void fixCenterPos(RDataFrame *df){
    int i = 0;
    for(int x = -3; x < 4; x+=1){
    for(int y = -3; y < 4; y+=1){
    for(int z = -3; z < 4; z+=1){

    
    auto *c = new TCanvas();
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);

    int bins = 300;
    double xmin = -1;
    double xmax = 1;


    auto d0 = df->Define("d0", "betaAlphaAngle0._betaAlphaAngle0");
    auto h0 = d0.Histo1D({"Stats", "data / efficiency", bins, xmin, xmax}, "d0");

    auto d1 = df->Define("d1", "betaAlphaAngle1._betaAlphaAngle1");
    auto h1 = d1.Histo1D({"Stats", "\\beta-\\alpha angle", bins, xmin, xmax}, "d1");

    
    h0->Add(&h1.GetValue()); 
    
    auto eff = new TH1F("stats", "Angle efficiency", bins, xmin, xmax);
    string filename = "try" + to_string(x) + to_string(y) + to_string(z);
    ifstream ifile("/home/anders/i257/build/effFiles/" + filename + ".csv");
    if (!ifile.is_open()) {
        std::cerr << "There was a problem opening the input file!\n";
        exit(1);//exit or do additional error checking
    }
    while (!ifile.eof())
    {
        double angle, weight;
        ifile >> angle >> weight;
        eff->Fill(angle, weight);
    }

    Double_t factor = 1.;
    eff->Scale(factor/eff->Integral(), "width");
    h0->Scale(factor/h0->Integral(), "width"); 

    auto xaxis = h0->GetXaxis();
    auto yaxis = h0->GetYaxis();
    xaxis->SetTitle("cos(a)");
    xaxis->CenterTitle();
    yaxis->SetTitle("count");
    yaxis->CenterTitle();
    eff->SetLineColor(kBlue);
    eff->SetLineWidth(3);
    h0->SetLineColor(kGreen);
    h0->SetLineWidth(3);
    
    h0->DrawClone("HIST"); 
    eff->DrawClone("HIST SAME");
    c->Modified();
    c->Update();
    auto hehehe = "/home/anders/i257/figures/centerCorrections/" + filename + ".png";
    auto save = hehehe.c_str();
    c->SaveAs(save);

    cout << save << endl;
    c->Close();
    ifile.close();
    i++;
    }
    }
    }
    

}

void momentumPlot(RDataFrame *df){
    auto c = new TCanvas();
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kViridis);
    int bins = 300;
    double xmin = 0;
    double xmax = 100000;

    auto d0 = df->Define("d0", "ptot._ptot");
    auto h0 = d0.Histo1D({"Stats", "data / efficiency", bins, xmin, xmax}, "d0");

    auto xaxis = h0->GetXaxis();
    auto yaxis = h0->GetYaxis();
    xaxis->SetTitle(" P_{total} [keV/c]");
    xaxis->CenterTitle();
    yaxis->SetTitle("count");
    yaxis->CenterTitle();
    h0->SetLineColor(kBlue);
    h0->SetLineWidth(3);


    h0->DrawClone();
    c->Modified();
    c->Update();
    c->SaveAs("/home/anders/i257/figures/ptotNoCut.pdf");
}

void singleAlphaSpectre(RDataFrame *df){
    auto c = new TCanvas();
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kViridis);
    auto data = df->Define("x", "E._E[0]");
    auto h = data.Histo1D({"Stats", "title", 300, 1, 8000}, "x");
    h->SetLineColor(kBlue);
    h->SetLineWidth(3);
    auto xaxis = h->GetXaxis();
    xaxis->SetTitle("E [keV]");
    xaxis->CenterTitle();
    h->DrawClone("Hist");
    c->SaveAs("/home/anders/i257/figures/singleAlpha.pdf");


}

void doubleAlphaSpectra(RDataFrame *df){
    auto c = new TCanvas();
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kViridis);
    auto data = df->Define("x", "E._E[0] + E._E[1]");
    auto data2 = df->Define("x", "E._E[1]");
    auto h = data.Histo1D({"Stats", "title", 300, 1, 18000}, "x");
    auto h2 = data2.Histo1D({"Stats", "title", 300, 1, 18000}, "x");

    h->SetLineColor(kBlue);
    h->SetLineWidth(3);
    h2->SetLineColor(kRed);
    h2->SetLineWidth(3);

    auto xaxis = h->GetXaxis();
    auto yaxis = h->GetYaxis();
    c->SetLogy();
    xaxis->SetTitle("E [keV]");
    xaxis->CenterTitle();
    
    h->DrawClone("Hist");
    // h2->DrawClone("HIST SAME");

    c->SaveAs("/home/anders/i257/figures/doubleAlpha.pdf");


}

void compareCuts(RDataFrame *dfAllCuts, RDataFrame *dfNoCuts, RDataFrame *dfAng, RDataFrame *dfMoment, RDataFrame *dfBetaMul, RDataFrame *dfAngAndMoment){
    auto c = new TCanvas();
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kViridis);

    auto all = dfAllCuts->Define("x", "E._E[0]");
    auto no = dfNoCuts->Define("x", "E._E[0]");
    auto ang = dfAng->Define("x", "E._E[0]");
    auto moment = dfMoment->Define("x", "E._E[0]");
    auto betaMul = dfBetaMul->Define("x", "E._E[0]");
    auto angAndMoment = dfAngAndMoment->Define("x", "E._E[0]");

    auto allH = all.Histo1D({"Stats", "titel", 300, 1, 8000}, "x");
    auto noH = no.Histo1D({"Stats", "titel", 300, 1, 8000}, "x");
    auto angH = ang.Histo1D({"Stats", "titel", 300, 1, 8000}, "x");
    auto momentH = moment.Histo1D({"Stats", "titel", 300, 1, 8000}, "x");
    auto betaMulH = betaMul.Histo1D({"Stats", "titel", 300, 1, 8000}, "x");
    auto angAndMomentH = angAndMoment.Histo1D({"Stats", "titel", 300, 1, 8000}, "x");


    allH->SetLineColor(kBlack);
    allH->SetLineWidth(3);
    noH->SetLineColor(kRed);
    noH->SetLineWidth(3);
    angH->SetLineColor(kGreen);
    angH->SetLineWidth(3);
    momentH->SetLineColor(kBlue);
    momentH->SetLineWidth(3);
    betaMulH->SetLineColor(kCyan);
    betaMulH->SetLineWidth(3);
    angAndMomentH->SetLineColor(kMagenta);
    angAndMomentH->SetLineWidth(3);

    auto xaxis = noH->GetXaxis();
    auto yaxis = noH->GetYaxis();
    yaxis->SetTitle("Counts");
    yaxis->CenterTitle();
    xaxis->SetTitle("E_{\\alpha 1} [keV]");
    xaxis->CenterTitle();
    yaxis->SetRange(0, 10);
    noH->DrawClone("HIST");
    allH->DrawClone("HIST SAME");
    angH->DrawClone("HIST SAME");
    momentH->DrawClone("HIST SAME");
    betaMulH->DrawClone("HIST SAME");
    angAndMomentH->DrawClone("HIST SAME");
    
    auto legend = new TLegend(0.55, 0.6, 0.88, 0.88);
    legend->AddEntry(noH->Clone(),"No cuts");
    legend->AddEntry(angH->Clone(),"Angular cut");
    legend->AddEntry(momentH->Clone(), "Momentum cut");
    legend->AddEntry(angAndMomentH->Clone(), "Angular and momentum cut");
    legend->AddEntry(betaMulH->Clone(), "Beta multiplicity cut");
    legend->AddEntry(allH->Clone(),"All cuts");
    legend->SetTextSize(0.03);
    legend->Draw();

    c->SetLogy();
    c->SaveAs("/home/anders/i257/figures/cutCompare.pdf");


}

void energyDifference(RDataFrame *df){
    auto c = new TCanvas();
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kViridis);
    auto data = df->Define("x", "E._E[0] - E._E[1]");
    auto h = data.Histo1D({"Stats", "Title", 300, -2000, 2000}, "x");
    h->SetLineColor(kBlue);
    h->SetLineWidth(3);
    h->DrawClone("hist");
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
    // TString filename = "/home/anders/i257/data/Li8/ONLY_2_ALPHAS_CUT_225_N102mlio.root";
    // TString filename = "/home/anders/i257/data/Li8/ONLY_ ANGULAR_CUT_225_N102mlio.root";
    // TString filename = "/home/anders/i257/data/Li8/ONLY_BETAMUL_CUT_225_N102mlio.root";
    // TString filename = "/home/anders/i257/data/Li8/ANG_AND_MOMENT_225_N102mlio.root";
    // TString filename = "/home/anders/i257/data/Li8/ONLY_EFF_CUT_225_N102mlio.root";
    //const char *input_file = argv[2];
    chain.Add(filename);
    RDataFrame df(chain);


    // Mulige konfigurationer:
    // No cut
    // Angular cut
    // Momentum cut (som en idiot åbenbart har kaldt eff??)
    // Beta Mul cut
    // Angular AND momentum cut
    // All cuts

    // No cut
    TChain CNoCut("tree");
    TString fileNoCut = "/home/anders/i257/data/Li8/ONLY_2_ALPHAS_CUT_225_N102mlio.root";
    CNoCut.Add(fileNoCut);
    RDataFrame dfNoCut(CNoCut);

    // Angular cut  
    TChain CAng("tree");
    TString fileAng = "/home/anders/i257/data/Li8/ONLY_ ANGULAR_CUT_225_N102mlio.root";
    CAng.Add(fileAng);
    RDataFrame dfAng(CAng);

    // Momentum cut (som en idiot åbenbart har kaldt eff??)
    TChain CMoment("tree");
    TString fileMoment = "/home/anders/i257/data/Li8/ONLY_EFF_CUT_225_N102mlio.root";
    CMoment.Add(fileMoment);
    RDataFrame dfMoment(CMoment);

    // Beta Mul cut
    TChain CBetaMul("tree");
    TString fileBetaMul = "/home/anders/i257/data/Li8/ONLY_BETAMUL_CUT_225_N102mlio.root";
    CBetaMul.Add(fileBetaMul);
    RDataFrame dfBetaMul(CBetaMul);

    // Angular AND momentum cut
    TChain CAngAndMoment("tree");
    TString fileAngAndMoment = "/home/anders/i257/data/Li8/ANG_AND_MOMENT_225_N102mlio.root";
    CAngAndMoment.Add(fileAngAndMoment);
    RDataFrame dfAngAndMoment(CAngAndMoment);

    // For legacy reasons hedder all cuts's dataframe df!




    // start a ROOT application window such that the plots can actually be shown
    TApplication *app = new TApplication("ROOT window", 0, 0);
    
    // Below here is all the plotting methods. Out comment the one you want to run. (Bad code dont judge!)

    // momentumPlot(&df);
    // fixCenterPos(&df);
    // singleDetectorEff(&df);
    // mexiHatTheory(&df);
    // mexiHatDetector(&df);
    // detectorEff(&df);
    // alphaEffeciency(&df);
    // cosang(&df);
    // EEfigure(&df);
    // betaAlphaDifferentEnergies(&df);
    // betaAlphaAngle(&df);
    // individualDetectorsBetaAlphaAngle(&df);
    // betaSpec(&df);
    // angEDiff(&df);
    // alphaAndBetaEnergy(&df);
    // singleAlphaSpectre(&df);
    // doubleAlphaSpectra(&df);
    // compareCuts(&df, &dfNoCut, &dfAng, &dfMoment, &dfBetaMul, &dfAngAndMoment);
    energyDifference(&df);
    app->Run(); // show all canvas
    return 0;
}
