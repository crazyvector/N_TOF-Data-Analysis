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
#include <sys/stat.h>   // mkdir
#include <sys/types.h>
#include <sstream>

//-----------------------------------------------------------
// Helper function to recursively create directories
//-----------------------------------------------------------
void createDirectories(const std::string &path) {
    std::stringstream ss(path);         // split the path by '/'
    std::string item;
    std::string currentPath = ".";      // start from current directory

    while (std::getline(ss, item, '/')) {  // iterate over each folder in path
        if (item.empty()) continue;       // skip empty strings (extra slashes)
        currentPath += "/" + item;        // append current folder
        mkdir(currentPath.c_str(), 0777); // create folder if it doesn't exist, 0777 = full permissions
    }
}

//-----------------------------------------------------------
// Function to convert a detector channel to energy
//-----------------------------------------------------------
bool EnergyCalibration(Float_t channel, Int_t Detector_ID,
                       const std::map<int, std::tuple<double, double, double>> &calib_params,
                       float &energy_out)
{
    static std::set<int> warned_detectors; // keeps track of detectors we already warned about

    // Find calibration parameters for this detector
    auto it = calib_params.find(Detector_ID);
    if (it == calib_params.end()) {
        // Warn once if detector calibration is missing
        if (warned_detectors.find(Detector_ID) == warned_detectors.end()) {
            std::cerr << "⚠️  Warning: Detector ID " << Detector_ID
                      << " not found in calibration parameters.\n";
            warned_detectors.insert(Detector_ID);
        }
        return false;  // calibration failed
    }

    // Extract calibration parameters (slope, intercept, quadratic term)
    double slope = std::get<0>(it->second);
    double intercept = std::get<1>(it->second);
    double quad = std::get<2>(it->second);

    // Apply quadratic calibration formula: E = intercept + slope*channel + quad*channel^2
    energy_out = intercept + slope * channel + quad * channel * channel;
    return true;
}

//-----------------------------------------------------------
// Load calibration parameters from a text file
//-----------------------------------------------------------
std::map<int, std::tuple<double, double, double>> LoadCalibrationParameters(
    const std::string &filename, const std::string &fit_type)
{
    std::map<int, std::tuple<double, double, double>> params;
    std::ifstream infile(filename);

    if (!infile.is_open()) {
        std::cerr << "❌ Error: Cannot open calibration parameters file: " << filename << std::endl;
        return params; // return empty map
    }

    int id;
    double slope, intercept, quad;
    std::string type;

    // Read file line by line: detectorID, slope, intercept, quad, type
    while (infile >> id >> slope >> intercept >> quad >> type) {
        if (type == fit_type) {
            // Keep only the first occurrence for this detector and fit type
            if (params.find(id) == params.end())
                params[id] = std::make_tuple(slope, intercept, quad);
        }
    }

    infile.close();
    return params;
}

//-----------------------------------------------------------
// Main function to create 2D histograms (Neutron Energy vs Gamma Energy)
//-----------------------------------------------------------
void plot2D(const std::vector<TString> &RootFileNrs, TString Detector,
            Int_t TotalDetector, TString Source, std::string fit_type,
            int T_MIN, int T_MAX, int ndec, int N_BPDEC, int nbinsX, double step,
            float ymin, float ymax, int n_bin_y, double binWidth,
            const std::string &calibrationFile, const std::string &rootFolder,
            double PsInt_threshold, double distance)
{
    // Create output folder for histograms
    std::string folder = "2D_Plots/";
    createDirectories(folder.c_str());

    //-------------------------------------------------------
    // Create a TChain to combine multiple ROOT files
    //-------------------------------------------------------
    TChain *chain = new TChain(Detector); // detector = tree name in ROOT files

    for (const auto &RootFileNr : RootFileNrs) {
        TString FName = rootFolder + "run" + RootFileNr + ".root";
        std::cout << "Adding file: " << FName << std::endl;

        // Check if file exists
        if (gSystem->AccessPathName(FName)) {
            std::cerr << "⚠️  Warning: File not found " << FName << std::endl;
            continue;
        }

        chain->Add(FName);  // add the file to TChain
    }

    //-------------------------------------------------------
    // Load calibration parameters from file
    //-------------------------------------------------------
    std::map<int, std::tuple<double, double, double>> calib_params =
        LoadCalibrationParameters(calibrationFile, fit_type);

    if (calib_params.empty()) {
        std::cerr << "❌ Error: Calibration parameters missing or invalid for fit type: "
                  << fit_type << std::endl;
        return;
    }

    //-------------------------------------------------------
    // Set variables to read from ROOT tree branches
    //-------------------------------------------------------
    Int_t detectorID;   // detector number
    Float_t channel;    // raw detector channel
    Float_t PsInt;      // pulse intensity
    Double_t TOF, TF;   // time-of-flight and time flash

    chain->SetBranchAddress("detn", &detectorID);        // detector ID branch
    chain->SetBranchAddress("amp", &channel);           // amplitude / channel branch
    chain->SetBranchAddress("tof", &TOF);               // TOF branch
    chain->SetBranchAddress("tflash", &TF);            // TF branch
    chain->SetBranchAddress("PulseIntensity", &PsInt); // pulse intensity

    //-------------------------------------------------------
    // Define X-axis bin edges (logarithmic)
    //-------------------------------------------------------
    std::vector<double> x_edges(nbinsX + 1);
    for (int i = 0; i <= nbinsX; i++)
        x_edges[i] = pow(10., T_MIN + i * step);

    Long64_t nEntries = chain->GetEntries(); // total number of events in TChain

    Long64_t progress_step = nEntries / 100; // progress every 1%
    if (progress_step == 0) progress_step = 1;

    //-------------------------------------------------------
    // Use a map to store histograms for each detector
    //-------------------------------------------------------
    std::map<int, TH2F *> hist;

    //-------------------------------------------------------
    // Loop over all events
    //-------------------------------------------------------
    for (Long64_t i_entry = 0; i_entry < nEntries; i_entry++) {
        if (i_entry % progress_step == 0) {
            std::cout << "Processing: " << (100 * i_entry / nEntries)
                      << "% (" << i_entry << "/" << nEntries << ")\n";
        }

        chain->GetEntry(i_entry); // load event

        //---------------------------------------------------
        // Skip events with pulse intensity below threshold
        //---------------------------------------------------
        if (PsInt <= PsInt_threshold) continue;

        //---------------------------------------------------
        // Convert channel to gamma energy using calibration
        //---------------------------------------------------
        float energy;
        if (!EnergyCalibration(channel, detectorID, calib_params, energy))
            continue; // skip if calibration fails

        //---------------------------------------------------
        // Calculate neutron energy from TOF
        //---------------------------------------------------
        double delta_t = TOF - TF;
        if (delta_t == 0) continue; // avoid division by zero

        float En = pow((72.2977 * distance / (delta_t * 1e-3)), 2); // neutron energy in eV
        float Eg = energy;                                          // gamma energy in keV

        //---------------------------------------------------
        // Create histogram for this detector if it doesn't exist
        //---------------------------------------------------
        if (hist.find(detectorID) == hist.end()) {
            TString hname = Form("histogram_%d", detectorID);
            hist[detectorID] = new TH2F(hname, "2D Histogram",
                                        nbinsX, x_edges.data(),
                                        n_bin_y, ymin, ymax);
        }

        //---------------------------------------------------
        // Fill histogram with (Neutron Energy, Gamma Energy)
        //---------------------------------------------------
        hist[detectorID]->Fill(En, Eg);
    }

    std::cout << "Processing: 100% (" << nEntries << "/" << nEntries << ")\n" << std::endl;

    //-------------------------------------------------------
    // Write histograms to output ROOT file
    //-------------------------------------------------------
    TString outfile_name = folder + "histograma_" + Detector + "_" + Source + ".root";
    TFile *outfile = new TFile(outfile_name, "RECREATE");

    for (auto &pair : hist) {
        TH2F *h = pair.second;
        if (h->GetEntries() > 0) {
            h->SetTitle(Form("Neutron Energy vs Gamma Energy - Detector %d", pair.first));
            h->GetXaxis()->SetTitle("Neutron Energy (eV)");
            h->GetYaxis()->SetTitle("Gamma Energy (KeV)");
            h->SetStats(0);           // disable statistics box
            h->SetContour(1000);      // smooth color contours

            outfile->cd();
            h->Write();               // save histogram

            //---------------------------------------------------
            // Optionally save PNG image for visualization
            //---------------------------------------------------
            TString imageName = Form("%shistograma_%s_%s_detector%d.png",
                                     folder.c_str(), Source.Data(), Detector.Data(), pair.first);
            TCanvas c;                 // local canvas
            gStyle->SetPalette(kRainBow); // color palette
            c.SetLogx();               // log scale X-axis
            h->Draw("COLZ");           // 2D color plot
            c.SaveAs(imageName);       // save PNG
        }

        delete h; // free memory
    }

    outfile->Close();
    std::cout << "✅ Histograms saved in file: " << outfile_name << std::endl;

    delete chain; // clean up TChain
}
