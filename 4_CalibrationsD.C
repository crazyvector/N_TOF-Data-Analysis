///////////////////////////////////////////////////
// Gabriel
// 01-08-2023 (modified for energy)
// Modified: 24-09-2025 - accepts parameters from JSON
///////////////////////////////////////////////////

#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TH1F.h>
#include <TChain.h>

#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <map>
#include <tuple>
#include <set>

Float_t EnergyCalibration(Float_t channel, Int_t Detector_ID, const std::map<int, std::tuple<double, double, double>> &calib_params)
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
        return 0.0;
    }

    double slope = std::get<0>(it->second);
    double intercept = std::get<1>(it->second);
    double quad = std::get<2>(it->second);

    return intercept + slope * channel + quad * channel * channel;
}

std::map<int, std::tuple<double, double, double>> LoadCalibrationParameters(TString filename, const std::string &fit_type)
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
            if (params.find(id) == params.end())
                params[id] = std::make_tuple(slope, intercept, quad);
        }
    }

    infile.close();
    return params;
}

void calibrations_energy(const std::vector<TString> &RootFileNrs, TString Detector, Int_t TotalDetector, TString Source, const std::string &fit_type, TString inputFolder, TString outputFolder, Int_t xminHist, Int_t xmaxHist, Int_t binWidth)
{
    // Loads calibration parameters
    auto calib_params = LoadCalibrationParameters(Form("%scalibration_parameters.txt", outputFolder.Data()), fit_type);
    if (calib_params.empty())
    {
        std::cerr << "❌ Error: No calibration parameters loaded.\n";
        return;
    }

    // Creates TChain
    TChain *chain = new TChain(Detector);
    for (const auto &RootFileNr : RootFileNrs)
    {
        TString FName = inputFolder + RootFileNr + ".root";
        chain->Add(FName);
    }

    if (chain->GetEntries() == 0)
    {
        std::cout << "❌ Error: No entries found in the chain.\n";
        delete chain;
        return;
    }

    // Variables
    Float_t channel = 0.0;
    Int_t detectorID;
    const Int_t xmin = xminHist;
    const Int_t xmax = xmaxHist;
    const Float_t bin_width = binWidth;
    const Float_t n_bin_x = (xmax - xmin) / bin_width;

    chain->SetBranchAddress("detn", &detectorID);
    chain->SetBranchAddress("amp", &channel);

    // Vector of histograms (in energy)
    std::vector<TH1F *> histo_energy(TotalDetector + 1);
    for (Int_t i = 1; i <= TotalDetector; i++)
    {
        TString name = Form("histograma_Energy_%d", i);
        histo_energy[i] = new TH1F(name, "Energy Spectrum", n_bin_x, xmin, xmax);
    }

    // Loop over all entries
    Long64_t nEntries = chain->GetEntries();
    for (Long64_t i_entry = 0; i_entry < nEntries; i_entry++)
    {
        chain->GetEntry(i_entry);
        if (detectorID >= 1 && detectorID <= TotalDetector)
        {
            Float_t energy = EnergyCalibration(channel, detectorID, calib_params);
            histo_energy[detectorID]->Fill(energy);
        }
    }

    // Scrie histogramele intr-un fisier ROOT
    TString outName;
    if (RootFileNrs.size() == 1)
        outName = outputFolder + "histograma_" + RootFileNrs[0] + "_" + Detector + "_" + Source + "_energy.root";
    else
        outName = outputFolder + "Rebin_" + Source + "_summed_histograms_energy.root";

    TFile *outfile = new TFile(outName, "RECREATE");
    for (Int_t i = 1; i <= TotalDetector; i++)
    {
        if (histo_energy[i])
            histo_energy[i]->Write();
    }
    outfile->Close();

    std::cout << "✅ Energy histograms saved in file: " << outName << "\n";

    delete chain;
}