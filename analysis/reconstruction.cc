#include "GetAngles.hh"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include <tuple>




int FillMultiEhist(TH1D* he, TH1D* hm, TH2D* h, const char* fname = "MNuball1.root") {
    std::map<int, std::map<std::vector<double>, double>> nuball = GetAngleDict(fname);
    int nofEvents = nuball.size();
    std::cout << nuball.size() << std::endl;
    for (const auto& n : nuball){
        double E = 0;
        int M = 0;
        for (const auto& a : n.second) {
            E += a.second;
            M++;
        }
        // std::cout << E << "   " << M << "\n";
        he->Fill(E);
        hm->Fill(M);
        h->Fill(E, M);
    }

    return nofEvents;
}

std::tuple<TH1D, TH1D, TH2D> GetHistGun(const char* fname = "MNuball.root") {
    TFile *inputf = new TFile(fname);
    double E;
    int M, evtID;
    TTree* t = (TTree*)inputf->Get("gun");
    t->SetBranchAddress("Energy", &E);
    t->SetBranchAddress("Multiplicity", &M);
    t->SetBranchAddress("EvtID", &evtID);
    TH1D *he = new TH1D("ghe", "Generated energy distribution;Energy, [MeV]",
                                                     40, 0., 40 /* 60, 0, 2 */);
    TH1D *hm = new TH1D("ghm", "Generated multiplicity distribution;Multiplicity", 
                                                        40, 0, 40 /* 15, 0, 15 */);
    TH2D *h = new TH2D("gh", "Generated Energy&Multiplicity distribution;Energy, [MeV];Multiplicity",
                                         40, 0., 40., 40, 0, 40/* 60, 0, 2, 15, 0, 15 */ );
    // std::map<int, int> m;
    for (int i = 0; i < t->GetEntries(); i++){
        t->GetEntry(i);
        he->Fill(E); 
        hm->Fill(M);
        h->Fill(E, M);
        // m[evtID]++;
    }

    // nofEvents = t->GetEntries();

    return {*he, *hm, *h};

}

void reconstruction(const char* fname = "gun.root"){

    const char* outfn = "reconstruction.root";
    TH1D *he = new TH1D("he", "Reconstructed energy distribution;Energy[MeV]",
                                                     40, 0., 40 /* 60, 0, 2 */);
    TH1D *hm = new TH1D("hm", "Reconstructed multiplicity distribution;Multiplicity", 
                                                        40, 0, 40 /* 15, 0, 15 */);
    TH2D *h = new TH2D("h", "Reconstructed Energy&Multiplicity distribution;Energy[MeV];Multiplicity",
                                         40, 0., 40., 40, 0, 40/* 60, 0, 2, 15, 0, 15 */ );
    int rN = FillMultiEhist(he, hm, h, fname);
    
    auto hists = GetHistGun(fname);
    TFile f(fname);
    // TTree* t = (TTree*)f.Get("gun");
    // int gN = t->GetEntries();
    // int w = gN - rN;
    // for (int i = 0; i < w; i++) hm->Fill(0);
    TFile *of = new TFile(outfn, "RECREATE");
    he->Write();
    hm->Write();
    h->Write();
    std::tuple<TH1D, TH1D, TH2D> rec_hists(*he, *hm, *h);
    std::get<0>(hists).Write();
    std::get<1>(hists).Write();
    std::get<2>(hists).Write();
    of->Close();

    TCanvas *c = new TCanvas("c", "c");

    c->SaveAs("reconstruction.pdf[");

    c->Divide(2, 1);
    c->cd(1);
    he->Draw();
    c->cd(2);
    std::get<0>(hists).Draw();
    c->SaveAs("reconstruction.pdf");
    c->Clear();

    c->Divide(2, 1);
    c->cd(1);
    hm->Draw();
    c->cd(2);
    std::get<1>(hists).Draw();
    c->SaveAs("reconstruction.pdf");
    c->Clear();

    c->Divide(2, 1);
    c->cd(1);
    h->Draw("colz");
    c->cd(2);
    std::get<2>(hists).Draw("colz");
    c->SaveAs("reconstruction.pdf");

    c->Clear();
    c->SaveAs("reconstruction.pdf]");
    
    
    
    
    
}