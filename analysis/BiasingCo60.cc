#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include <set>
#include <assert.h>
#include <iostream>
#include <algorithm>
#include "GetAngles.hh"


#include "TH1D.h"
#include "TCanvas.h"

#include "TApplication.h"

const char* fname = "Co60.root";
// const char* fname = "Nuballt.root";
double deg = M_PI/180;

bool IsPhotoPeak (double E, double eps = 0.01) {
    double E1 = 1.1732;
    double E2 = 1.3325;
    return (abs(E1 - E) < eps) || (abs(E2 - E) < eps);
}

void BiasingCo60 (const char* name = fname) {

    double Eg1 = 1.1732;
    double Eg2 = 1.3325;
    TFile *f = new TFile(name);
    int id_1, id_2;
    double x1, x2, y1, y2, z1, z2; 
    TTree *p1 = (TTree*)f->Get("p1");
    p1->SetBranchAddress("x", &x1);
    p1->SetBranchAddress("y", &y1);
    p1->SetBranchAddress("z", &z1);
    p1->SetBranchAddress("EvtID", &id_1);
    TTree *p2 = (TTree*)f->Get("p2");
    p2->SetBranchAddress("x", &x2);
    p2->SetBranchAddress("y", &y2);
    p2->SetBranchAddress("z", &z2);
    p2->SetBranchAddress("EvtID", &id_2);
    double *alpha0 = new double[p1->GetEntries()];
    std::map<int, TVector3> v1;
    std::map<int, TVector3> v2;
    std::cout << "Entries: " << p1->GetEntries() << std::endl;
    for (int i = 0; i < p1->GetEntries(); ++i) {
        p1->GetEntry(i);
        p2->GetEntry(i);
        v1[id_1] = TVector3(x1, y1, z1);
        v2[id_2] = TVector3(x2, y2, z2);   
    }
    TH1D *h0 = new TH1D("h0", "generated distribution", 18, 0., 180.);
    std::cout << "starting filling histograms\n";
    for (int i = 0; i < v1.size(); ++i) {
        if (v1.find(i) == v1.end() || v2.find(i) == v2.end()) {
            std::cout << "something occurs\n";
        }
        alpha0[i] = v1[i].Angle(v2[i]);
        h0->Fill(alpha0[i]/deg);
    }
    assert(v1.size() == v2.size());

    std::cout << "first histo was filled\n";

    TCanvas *c = new TCanvas();
    c->Divide(2,1);
    c->cd(1);
    h0->Draw();

    std::map<int, std::map<std::vector<double>, double>> nuball = GetbAngleDict(fname);
    TH1D *h = new TH1D("h", "counts", 18, 0., 180.);
    

    for (const auto & e:nuball) {
        int EvtID = e.first;
        int M = 0;
        for (const auto& el: nuball[EvtID]) {
            double Energy_per_det = el.second;
            if (IsPhotoPeak(Energy_per_det)) M++;
        }
        if (M == 2) h->Fill(alpha0[EvtID]/deg);
        if (M > 2) std::cout 
                << "Warning: Multiplicity > 2 !\n";
        

    }
    std::cout << "the second histo was filled\n";
    delete [] alpha0;
    // double *ErrorSet = new double[nuball.size()];
    int nBins = h->GetNbinsX();
    int nBins0 = h0->GetNbinsX();
    assert(nBins0 == nBins);
    TH1D *hp = new TH1D("hp", "probability to detect both gammas", 18, 0., 180.);
    std::cout << "Detections:\n";
    for (int i = 1; i <= nBins; ++i) {


        double Y = h->GetBinContent(i)/h0->GetBinContent(i);
        double deltaX = h->GetBinError(i);
        if (h->GetBinContent(i) == 0) {deltaX = (int)0; Y = (int)0;}
        assert(h->GetBinError(i) == sqrt(h->GetBinContent(i)));
        double sigmaX = deltaX/h->GetBinContent(i);
        if (h->GetBinContent(i) == 0) sigmaX = 0;
        double sigmaY = Y * sigmaX;
        hp->Fill(h->GetBinCenter(i), Y);
        hp->SetBinError(i, sigmaY);
        std::cout << h->GetBinCenter(i) << "\t" << "Genereted = " << h0->GetBinContent(i) << "\t"
                  << "X = " << h->GetBinContent(i) << "\t" 
                  << "deltaX " << deltaX << "\t" << "\t"
                  << "precision = " << sigmaX << "\t"
                  << "Y = " << Y << "\t" << "Er = " << sigmaY << "\n";
    }
    std::cout << "third histo were calculated\n";

    c->cd(2);
    h->Draw("E");
    
    
    // *hp = *h/(*h0);
    // std::cout << "dividation completed\n";
    // ErrorSet[0] = 0;
    // hp->SetTitle("probability to detect both gammas");
    TCanvas *nc = new TCanvas();
    // hp->SetError(ErrorSet);
    // delete [] ErrorSet;
    // hp->SetBinError(0, 0);
    hp->GetYaxis()->SetRangeUser(0, 0.007);
    hp->Draw("E");
    c->SaveAs("Co60counts0.png");
    
    nc->SaveAs("Co60probability0.png");

    for (int i = 1; i <= hp->GetNbinsX(); ++i) {
        std::cout << hp->GetBinCenter(i) << "\t" << hp->GetBinContent(i) << "\n";
    }
    
    
    
        

    
}

// int main (int argc, char **argv) {
//     TApplication app("app", &argc, argv);
//     // if (argc > 1) {
//     //     BiasingCo60(argv[1]);
//     // } else 
//     // {BiasingCo60();}
//     BiasingCo60("Nuball1.root");
//     app.Run();
//     return 0;
// }