#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <set>

#include "TTree.h"
#include "TFile.h"




std::map<int, std::vector<double>> GetAngles()
{

    std::map<int, std::vector<double>> angles;

    double AngleBoardBGOmin, AngleBoardBGOMax, AngleBoardmin, AngleBoardmax;
    std::set<std::string> file_names =
        {"PhaseIPosition.dat", "PhaseIBGOPosition.dat", "CloverBGOPosition.dat", "CloverPosition.dat", "LaBr3Position.dat"};
    std::ifstream input;
    int line_ind = 1;
    for (const auto &name : file_names)
    {
        input = std::ifstream(name);
        if (input)
        {
            std::string line;
            while (getline(input, line))
            {
                // std::cout << line << std::endl;
                std::stringstream dat(line);
                double rad, theta, phi;
                int label;
                if (line_ind == 3)
                {
                    //    std::string s;
                    //    dat >> s;
                    dat >> rad >> theta >> phi;
                    //    std::cout << line << "\t" << rad << " " << theta << " " << phi << "\n";
                    getline(input, line);
                    dat = std::stringstream(line);
                    ++line_ind;
                }
                if (line_ind == 4)
                {
                    dat >> label;
                    angles[label] = {rad, theta, phi};
                    line_ind = 0;
                    //    std:: cout << "label: " << label << "coord: "  << rad << " " << phi << " " << theta << std ::endl;
                }

                ++line_ind;
            }
        }
        else
        {
            std::cout << "<<< error in " << name << std::endl;
        }
    }

    for (int l = 24; l <= 162; l += 6)
    {

        angles[l - 1] = angles[l];
        auto sphc = angles[l + 1];
        for (int i = l + 2; i < l + 5; i++)
        {
            angles[i] = sphc;
        }
        if (angles[l] != angles[l + 1])
        {
            std::cout << "Error!!! angles[l] != angles[l+1]" << std::endl;
            std::cout << "l = " << l << std::endl
                      << std::endl;
        }
    }

    std::ofstream output("angels.dat");
    for (const auto &e : angles)
    {
        output << e.first << "\t" << e.second[0] << " " << e.second[1] << " " << e.second[2] << "\n";
    }

    return angles;
}

std::map<int, std::map<std::vector<double>, double>> GetAngleDict(const char* fname = "MNuball.root") {
    std::map<int, std::map<std::vector<double>, double>> res;
    TFile* f = new TFile(fname);
    TTree* t = (TTree*)f->Get("Nuball");
    int EvtID, l;
    double E;
    t->SetBranchAddress("Energy", &E);
    t->SetBranchAddress("Label", &l);
    t->SetBranchAddress("EventNumber", &EvtID);
    std::cout << t->GetEntries() << std::endl;
    auto angles = GetAngles();
    for (int i = 0; i < t->GetEntries(); ++i){
        t->GetEntry(i);
        res[EvtID][angles[l]] = res[EvtID][angles[l]] + E;
    }
    delete f;
    return res;
}
std::map<int, std::map<std::vector<double>, double>> GetbAngleDict(const char* fname) {
    std::map<int, std::map<std::vector<double>, double>> res;
    TFile* f = new TFile(fname);
    TTree* t = (TTree*)f->Get("Nuball");
    int EvtID, l, type;
    double E;
    t->SetBranchAddress("Energy", &E);
    t->SetBranchAddress("Label", &l);
    t->SetBranchAddress("EventNumber", &EvtID);
    t->SetBranchAddress("Type", &type);
    std::cout << t->GetEntries() << std::endl;
    auto angles = GetAngles();
    for (int i = 0; i < t->GetEntries(); ++i){
        t->GetEntry(i);
        if (type == 1)
        res[EvtID][angles[l]] = res[EvtID][angles[l]] + E;
    }
    return res;
}
std::map<int, std::map<int, bool>> GetTypeDict(const char* fname = "MNuball.root", int itype = 1) {
    std::map<int, std::map<int, bool>> res;

    TFile* f = new TFile(fname);
    TTree* t = (TTree*)f->Get("Nuball");
    int EvtID, l, type;
    double E;
    t->SetBranchAddress("Energy", &E);
    t->SetBranchAddress("Label", &l);
    t->SetBranchAddress("EventNumber", &EvtID);
    t->SetBranchAddress("Type", &type);

    // std::cout << t->GetEntries() << std::endl;

    
    for (int i = 0; i < t->GetEntries(); ++i){
        t->GetEntry(i);
        res[EvtID][l] = (type == itype);
    }
    return res;
}


std::map<int, int> GetMDict(const char* fname, int itype = 1, double gunEnergy = 1.) {
    std::map<int, int> nsuccess;
    TFile* f = new TFile(fname);
    TTree* t = (TTree*)f->Get("Nuball");
    int EvtID, l, type;
    double E;
    t->SetBranchAddress("Energy", &E);
    // t->SetBranchAddress("Label", &l);
    t->SetBranchAddress("EventNumber", &EvtID);
    t->SetBranchAddress("Type", &type);
    
    for (int i = 0; i < t->GetEntries(); ++i){
        t->GetEntry(i);
        nsuccess[EvtID] += ((abs(E - gunEnergy) < 0.01) && (type == itype));
    }
    return nsuccess;
}

std::map<int, std::map<int, double>> GetLabelDict(const char* fname = "Nuball.root") {
    std::map<int, std::map<int, double>> res;
    TFile* f = new TFile(fname);
    TTree* t = (TTree*)f->Get("Nuball");
    int EvtID, l;
    double E;
    t->SetBranchAddress("Energy", &E);
    t->SetBranchAddress("Label", &l);
    t->SetBranchAddress("EventNumber", &EvtID);
    for (int i = 0; i < t->GetEntries(); ++i){
        t->GetEntry(i);
        res[EvtID][l] = E;
    }
    return res;
}