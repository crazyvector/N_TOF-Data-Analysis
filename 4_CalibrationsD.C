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

// ===========================================================
// Function: EnergyCalibration
// Calculates energy from channel for a given detector
// using calibration parameters (slope, intercept, quadratic term)
// ===========================================================
Float_t EnergyCalibration(Float_t channel, Int_t Detector_ID, const std::map<int, std::tuple<double, double, double>> &calib_params)
{
    static std::set<int> warned_detectors; // keep track of detectors missing calibration

    auto it = calib_params.find(Detector_ID);
    if (it == calib_params.end())
    {
        // Only warn once per missing detector
        if (warned_detectors.find(Detector_ID) == warned_detectors.end())
        {
            std::cerr << "⚠️  Warning: Detector ID " << Detector_ID << " not found in calibration parameters.\n";
            warned_detectors.insert(Detector_ID);
        }
        return 0.0; // return zero if calibration parameters are missing
    }

    // Extract slope, intercept, and quadratic term
    double slope = std::get<0>(it->second);
    double intercept = std::get<1>(it->second);
    double quad = std::get<2>(it->second);

    // Return calibrated energy
    return intercept + slope * channel + quad * channel * channel;
}

// ===========================================================
// Function: LoadCalibrationParameters
// Reads calibration parameters from a text file and returns a map
// Only parameters matching the requested fit_type are loaded
// ===========================================================
std::map<int, std::tuple<double, double, double>> LoadCalibrationParameters(TString filename, const std::string &fit_type)
{
    std::map<int, std::tuple<double, double, double>> params;
    std::ifstream infile(filename);

    if (!infile.is_open())
    {
        std::cerr << "❌ Error: Cannot open calibration parameters file: " << filename << std::endl;
        return params; // return empty map if file cannot be opened
    }

    int id;
    double slope, intercept, quad;
    std::string type;

    while (infile >> id >> slope >> intercept >> quad >> type)
    {
        // Only store parameters matching the requested fit type
        if (type == fit_type)
        {
            if (params.find(id) == params.end())
                params[id] = std::make_tuple(slope, intercept, quad);
        }
    }

    infile.close();
    return params;
}

// ===========================================================
// Function: calibrations_energy
// Converts channel histograms to energy histograms using calibration parameters
// Supports multiple ROOT files and detectors
// ===========================================================
void calibrations_energy(const std::vector<TString> &RootFileNrs, TString Detector, Int_t TotalDetector, TString Source, const std::string &fit_type, TString inputFolder, TString outputFolder, Int_t xminHist, Int_t xmaxHist, Int_t binWidth)
{
    // Load calibration parameters from file
    auto calib_params = LoadCalibrationParameters(Form("%scalibration_parameters.txt", outputFolder.Data()), fit_type);
    if (calib_params.empty())
    {
        std::cerr << "❌ Error: No calibration parameters loaded.\n";
        return; // exit if no calibration parameters
    }

    // Create a TChain to combine multiple ROOT files
    TChain *chain = new TChain(Detector);
    for (const auto &RootFileNr : RootFileNrs)
    {
        TString FName = inputFolder + RootFileNr + ".root";
        chain->Add(FName); // add each ROOT file to the chain
    }

    if (chain->GetEntries() == 0)
    {
        std::cout << "❌ Error: No entries found in the chain.\n";
        delete chain;
        return; // exit if chain is empty
    }

    // Define variables for branches
    Float_t channel = 0.0;
    Int_t detectorID;

    // Histogram parameters
    const Int_t xmin = xminHist;
    const Int_t xmax = xmaxHist;
    const Float_t bin_width = binWidth;
    const Float_t n_bin_x = (xmax - xmin) / bin_width;

    // Set branch addresses for reading data
    chain->SetBranchAddress("detn", &detectorID); // detector ID
    chain->SetBranchAddress("amp", &channel);     // signal amplitude (channel)

    // Create a vector of TH1F histograms for each detector
    std::vector<TH1F *> histo_energy(TotalDetector + 1);
    for (Int_t i = 1; i <= TotalDetector; i++)
    {
        TString name = Form("histograma_Energy_%d", i);
        histo_energy[i] = new TH1F(name, "Energy Spectrum", n_bin_x, xmin, xmax);
    }

    // Loop over all entries in the TChain
    Long64_t nEntries = chain->GetEntries();
    for (Long64_t i_entry = 0; i_entry < nEntries; i_entry++)
    {
        chain->GetEntry(i_entry);
        if (detectorID >= 1 && detectorID <= TotalDetector)
        {
            // Convert channel to energy using calibration parameters
            Float_t energy = EnergyCalibration(channel, detectorID, calib_params);
            histo_energy[detectorID]->Fill(energy); // fill histogram
        }
    }

    // Determine output ROOT filename
    TString outName;
    if (RootFileNrs.size() == 1)
        outName = outputFolder + "histograma_" + RootFileNrs[0] + "_" + Detector + "_" + Source + "_energy.root";
    else
        outName = outputFolder + "Rebin_" + Source + "_summed_histograms_energy.root";

    // Create output ROOT file and save histograms
    TFile *outfile = new TFile(outName, "RECREATE");
    for (Int_t i = 1; i <= TotalDetector; i++)
    {
        if (histo_energy[i])
            histo_energy[i]->Write();
    }
    outfile->Close();

    std::cout << "✅ Energy histograms saved in file: " << outName << "\n";

    // Clean up
    delete chain;
}
