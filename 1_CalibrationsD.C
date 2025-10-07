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

///////////////////////////////////////////////////
// Recursive helper function for directory creation
///////////////////////////////////////////////////
void createDirectories(const std::string &path) {
    // Create a stringstream from the given path
    std::stringstream ss(path);
    std::string item;
    std::string currentPath = "."; // Start from the current directory (".")

    // Split the path by '/' and create directories step by step
    while (std::getline(ss, item, '/')) {
        if (item.empty()) continue; // skip empty parts (e.g., double slashes)
        currentPath += "/" + item;
        // mkdir() tries to create the directory; if it already exists, it does nothing
        mkdir(currentPath.c_str(), 0777); 
    }
}

///////////////////////////////////////////////////
// Main calibration function
// Builds histograms for each detector based on amplitude data
///////////////////////////////////////////////////
void calibrations(const std::vector<TString> &RootFileNrs, // list of ROOT file identifiers (filenames without extension)
                  TString Detector,                        // name of the detector branch in the tree
                  Int_t TotalDetector,                     // total number of detectors
                  TString Source,                          // name of the calibration source (used in output filename)
                  TString inputFolder,                     // folder containing input ROOT files
                  TString outputFolder,                    // folder where histograms will be saved
                  Int_t rebinFactor,                       // factor for histogram rebinning
                  Int_t xminHist,                          // lower x-axis limit for histograms
                  Int_t xmaxHist,                          // upper x-axis limit for histograms
                  Int_t binWidth)                          // bin width for histograms
{
    // Ensure that the output directory exists before saving histograms
    createDirectories(outputFolder.Data());

    ///////////////////////////////////////////////////
    // Create a TChain that merges all ROOT files together
    ///////////////////////////////////////////////////
    // Using TChain instead of TFile/TTree allows reading multiple runs easily
    TChain *chain = new TChain(Detector);

    // Loop through each file number provided in the list
    for (const auto &RootFileNr : RootFileNrs)
    {
        // Build the full path of each ROOT file
        TString FName = inputFolder + RootFileNr + ".root";
        // Add the file to the chain; ROOT will internally link all TTrees with the same name
        chain->Add(FName);
    }

    // Check if there are any entries to process; exit if not
    if (chain->GetEntries() == 0)
    {
        std::cout << "❌ Error: No entries found in the chain.\n";
        delete chain;
        return;
    }

    ///////////////////////////////////////////////////
    // Define temporary variables for reading data
    ///////////////////////////////////////////////////
    Float_t amplitude = 0.0; // temporary amplitude value for each entry
    Int_t detectorID;        // detector identifier (e.g., detn)
    Float_t channel;         // channel or amplitude value (raw data from tree)

    // Define histogram boundaries and binning
    const Int_t xmin = xminHist;
    const Int_t xmax = xmaxHist;
    const Float_t bin_width = binWidth;
    const Float_t n_bin_x = (xmax - xmin) / bin_width; // number of bins on x-axis

    ///////////////////////////////////////////////////
    // Connect tree branches to local variables
    ///////////////////////////////////////////////////
    chain->SetBranchAddress("detn", &detectorID); // detector ID
    chain->SetBranchAddress("amp", &channel);     // amplitude/channel value

    ///////////////////////////////////////////////////
    // Create one histogram for each detector
    ///////////////////////////////////////////////////
    std::vector<TH1F *> histo_channels(TotalDetector + 1); // index starts from 1 (not 0)

    for (Int_t i = 1; i <= TotalDetector; i++)
    {
        // Unique histogram name per detector
        TString name = Form("histograma_Amplitude_%d", i);
        // Create a histogram with uniform binning
        histo_channels[i] = new TH1F(name, "Amplitude", n_bin_x, xmin, xmax);
    }

    ///////////////////////////////////////////////////
    // Loop over all entries in all runs
    ///////////////////////////////////////////////////
    Long64_t nEntries = chain->GetEntries(); // total number of events
    for (Long64_t i_entry = 0; i_entry < nEntries; i_entry++)
    {
        // Load one entry from the chain into the local variables
        chain->GetEntry(i_entry);

        // Only fill histograms for valid detector IDs
        if (detectorID >= 1 && detectorID <= TotalDetector)
        {
            amplitude = static_cast<Float_t>(channel); // cast channel to float
            histo_channels[detectorID]->Fill(amplitude); // fill the corresponding histogram
        }
    }

    ///////////////////////////////////////////////////
    // Prepare output file name
    ///////////////////////////////////////////////////
    TString outName;
    if (RootFileNrs.size() == 1)
    {
        // Single file mode → keep specific run and detector info
        outName = outputFolder + "histograma_" + RootFileNrs[0] + "_" + Detector + "_" + Source + ".root";
    }
    else
    {
        // Multi-file mode → combined histograms
        outName = outputFolder + "Rebin_" + Source + "_summed_histograms.root";
    }

    ///////////////////////////////////////////////////
    // Save all histograms into a single ROOT file
    ///////////////////////////////////////////////////
    TFile *outfile = new TFile(outName, "RECREATE"); // overwrite if file exists

    for (Int_t i = 1; i <= TotalDetector; i++)
    {
        if (histo_channels[i])
        {
            // Apply rebinning only if requested (factor > 1)
            if (rebinFactor > 1) {
                histo_channels[i]->Rebin(rebinFactor);
            }

            // Write histogram to the output file
            histo_channels[i]->Write();
        }
    }

    // Close output file to save changes
    outfile->Close();

    ///////////////////////////////////////////////////
    // Final confirmation message
    ///////////////////////////////////////////////////
    std::cout << "✅ Histograms saved in file: " << outName << "\n";

    // Free memory (important when processing many runs)
    delete chain;
}
