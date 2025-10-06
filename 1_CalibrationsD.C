///////////////////////////////////////////////////
// Gabriel
// 01-08-2023
// Modificat: 24-09-2025 - accepta parametri din JSON
///////////////////////////////////////////////////

#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TChain.h>

#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

#include <sys/stat.h>
#include <sys/types.h>
#include <sstream>

// helper recursiv
void createDirectories(const std::string &path) {
    std::stringstream ss(path);
    std::string item;
    std::string currentPath = "."; // pornim din directorul curent

    while (std::getline(ss, item, '/')) {
        if (item.empty()) continue;
        currentPath += "/" + item;
        mkdir(currentPath.c_str(), 0777); // daca exista deja, nu da eroare
    }
}

void calibrations(const std::vector<TString> &RootFileNrs, TString Detector, Int_t TotalDetector, TString Source, TString inputFolder, TString outputFolder, Int_t rebinFactor, Int_t xminHist, Int_t xmaxHist, Int_t binWidth)
{
    createDirectories(outputFolder.Data());
    //mkdir(outputFolder.Data(), 0777);

    // Creeaza TChain in loc de TFile/TTree
    TChain *chain = new TChain(Detector);
    for (const auto &RootFileNr : RootFileNrs)
    {
        // Path-ul catre fisier este acum un parametru
        TString FName = inputFolder + RootFileNr + ".root";
        chain->Add(FName);
    }

    if (chain->GetEntries() == 0)
    {
        std::cout << "❌ Error: No entries found in the chain.\n";
        delete chain;
        return;
    }

    // Variabile temporare
    Float_t amplitude = 0.0;
    Int_t detectorID;
    Float_t channel;
    const Int_t xmin = xminHist;
    const Int_t xmax = xmaxHist;
    const Float_t bin_width = binWidth;
    const Float_t n_bin_x = (xmax - xmin) / bin_width;

    // Seteaza branch-urile
    chain->SetBranchAddress("detn", &detectorID);
    chain->SetBranchAddress("amp", &channel);

    // Creeaza vectorul de histograme
    std::vector<TH1F *> histo_channels(TotalDetector + 1);
    for (Int_t i = 1; i <= TotalDetector; i++)
    {
        TString name = Form("histograma_Amplitude_%d", i);
        histo_channels[i] = new TH1F(name, "Amplitude", n_bin_x, xmin, xmax);
    }

    // Loop peste toate intrarile din toate run-urile
    Long64_t nEntries = chain->GetEntries();
    for (Long64_t i_entry = 0; i_entry < nEntries; i_entry++)
    {
        chain->GetEntry(i_entry);
        if (detectorID >= 1 && detectorID <= TotalDetector)
        {
            amplitude = static_cast<Float_t>(channel);
            histo_channels[detectorID]->Fill(amplitude);
        }
    }

    // Scrie histogramele intr-un singur fisier ROOT
    TString outName;
    if (RootFileNrs.size() == 1)
    {
        outName = outputFolder + "histograma_" + RootFileNrs[0] + "_" + Detector + "_" + Source + ".root";
    }
    else
    {
        outName = outputFolder + "Rebin_" + Source + "_summed_histograms.root";
    }

    TFile *outfile = new TFile(outName, "RECREATE");

    for (Int_t i = 1; i <= TotalDetector; i++)
    {
        if (histo_channels[i])
        {
            // Valoarea de rebinning este acum un parametru
            if (rebinFactor > 1) {
                histo_channels[i]->Rebin(rebinFactor);
            }
            histo_channels[i]->Write();
        }
    }
    outfile->Close();

    std::cout << "✅ Histograms saved in file: " << outName << "\n";

    delete chain;
}
