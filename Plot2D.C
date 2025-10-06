#include <TString.h>
#include <TChain.h>
#include <TH2F.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TROOT.h>
#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <tuple>
#include <vector>
#include <string>
#include <cmath>

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

bool EnergyCalibration(Float_t channel, Int_t Detector_ID, const std::map<int, std::tuple<double, double, double>> &calib_params, float &energy_out)
{
    static std::set<int> warned_detectors;

    auto it = calib_params.find(Detector_ID);
    if (it == calib_params.end())
    {
        if (warned_detectors.find(Detector_ID) == warned_detectors.end())
        {
            std::cerr << "⚠️  Warning: Detector ID " << Detector_ID << " not found in calibration parameters.\n";
            warned_detectors.insert(Detector_ID);
        }
        return false;
    }

    double slope = std::get<0>(it->second);
    double intercept = std::get<1>(it->second);
    double quad = std::get<2>(it->second);

    energy_out = intercept + slope * channel + quad * channel * channel;
    return true;
}

std::map<int, std::tuple<double, double, double>> LoadCalibrationParameters(const std::string &filename, const std::string &fit_type)
{
    std::map<int, std::tuple<double, double, double>> params;
    std::ifstream infile(filename);

    if (!infile.is_open())
    {
        std::cerr << "❌ Error: Cannot open calibration parameters file: " << filename << std::endl;
        return params;
    }

    int id;
    double slope, intercept, quad;
    std::string type;

    while (infile >> id >> slope >> intercept >> quad >> type)
    {
        if (type == fit_type)
        {
            // Retine doar prima aparitie pentru acel detector cu acel tip de fit
            if (params.find(id) == params.end())
                params[id] = std::make_tuple(slope, intercept, quad);
        }
    }

    infile.close();
    return params;
}

void plot2D(const std::vector<TString> &RootFileNrs, TString Detector, Int_t TotalDetector, TString Source, std::string fit_type, int T_MIN, int T_MAX, int ndec, int N_BPDEC, int nbinsX, double step,
            float ymin, float ymax, int n_bin_y, double binWidth, const std::string &calibrationFile, const std::string &rootFolder, double PsInt_threshold, double distance)

{
    std::string folder = "2D_Plots/";
    createDirectories(folder.c_str());
    // if (gSystem->AccessPathName(folder.c_str()))
    // {
    //     gSystem->mkdir(folder.c_str(), true); // true = recursive
    // }

    // Creeaza TChain
    TChain *chain = new TChain(Detector);
    for (const auto &RootFileNr : RootFileNrs)
    {
        TString FName = rootFolder + "run" + RootFileNr + ".root";
        cout << "Adding file: " << FName << endl;
        if (gSystem->AccessPathName(FName))
        {
            std::cerr << "⚠️  Warning: Nu am găsit fișierul " << FName << std::endl;
            continue;
        }
        chain->Add(FName);
    }

    // Parametri de calibrare
    std::map<int, std::tuple<double, double, double>> calib_params = LoadCalibrationParameters(calibrationFile, fit_type);
    if (calib_params.empty())
    {
        std::cerr << "❌ Error: Calibration parameters missing or invalid for fit type: " << fit_type << std::endl;
        return;
    }

    // Variabile temporare
    Int_t detectorID;
    Float_t channel, PsInt;
    Double_t TOF, TF;

    // Setează branch-urile
    chain->SetBranchAddress("detn", &detectorID);
    chain->SetBranchAddress("amp", &channel);
    chain->SetBranchAddress("tof", &TOF);
    chain->SetBranchAddress("tflash", &TF);
    chain->SetBranchAddress("PulseIntensity", &PsInt);

    std::vector<double> x_edges(nbinsX + 1);
    for (int i = 0; i <= nbinsX; i++)
        x_edges[i] = pow(10., T_MIN + i * step);

    Long64_t nEntries = chain->GetEntries();

    Long64_t progress_step = nEntries / 100;
    if (progress_step == 0)
        progress_step = 1;

    // ---- OPTIMIZARE: folosim map pentru histograme ----
    std::map<int, TH2F *> hist;

    for (Long64_t i_entry = 0; i_entry < nEntries; i_entry++)
    {
        if (i_entry % progress_step == 0)
        {
            std::cout << "Processing: " << (100 * i_entry / nEntries) << "% (" << i_entry << "/" << nEntries << ")\n";
        }

        chain->GetEntry(i_entry);

        if (PsInt <= PsInt_threshold)
            continue; // skip evenimente cu intensitate mai mică de threshold

        float energy;
        if (!EnergyCalibration(channel, detectorID, calib_params, energy))
            continue;

        double delta_t = TOF - TF;
        if (delta_t == 0)
            continue; // skip ca să evită div by zero

        float En = pow((72.2977 * distance / (delta_t * 1e-3)), 2); // neutron energy in eV
        float Eg = energy;                                          // gamma energy in keV

        // Creeaza histograma doar daca nu exista deja
        if (hist.find(detectorID) == hist.end())
        {
            TString hname = Form("histogram_%d", detectorID);
            hist[detectorID] = new TH2F(hname, "2D Histogram",
                                        nbinsX, x_edges.data(),
                                        n_bin_y, ymin, ymax);
        }

        hist[detectorID]->Fill(En, Eg);
    }
    std::cout << "Processing: 100% (" << nEntries << "/" << nEntries << ")\n"
              << std::endl;

    // Fisier ROOT pentru output
    TString outfile_name = folder + "histograma_" + Detector + "_" + Source + ".root";
    TFile *outfile = new TFile(outfile_name, "RECREATE");

    // ---- Scrie fiecare histograma ----
    for (auto &pair : hist)
    {
        TH2F *h = pair.second;
        if (h->GetEntries() > 0)
        {
            h->SetTitle(Form("Neutron Energy vs Gamma Energy - Detector %d", pair.first));
            h->GetXaxis()->SetTitle("Neutron Energy (eV)");
            h->GetYaxis()->SetTitle("Gamma Energy (KeV)");
            h->SetStats(0);
            h->SetContour(1000);

            outfile->cd();
            h->Write();

            // Optional: genereaza imagine PNG doar dupa ce totul e scris in ROOT
            TString imageName = Form("%shistograma_%s_%s_detector%d.png", folder.c_str(), Source.Data(), Detector.Data(), pair.first);
            TCanvas c; // canvas local
            gStyle->SetPalette(kRainBow);
            c.SetLogx();
            h->Draw("COLZ");
            c.SaveAs(imageName);
        }

        delete h; // elibereaza RAM
    }

    outfile->Close();
    std::cout << "✅ Histograms saved in file: " << outfile_name << std::endl;

    delete chain;
}