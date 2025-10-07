#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TFitResult.h>
#include <TPaveText.h>
#include <TSpectrum.h>
#include <TF1.h>
#include <TMath.h>
#include <TLine.h>
#include <TText.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <utility>
#include <algorithm>
#include <iomanip>
#include <string>
#include <sstream>

///////////////////////////////////////////////////////////
// A simple structure to store information about a found peak
///////////////////////////////////////////////////////////
struct PeakInfo
{
    double position;   // X position (channel) of the peak
    double fwhm;       // Full Width at Half Maximum estimate
    double amplitude;  // Peak amplitude (height)
};

///////////////////////////////////////////////////////////
// Helper function to save a ROOT canvas into a multi-page PDF
///////////////////////////////////////////////////////////
void SaveCanvasToPDF(TCanvas *canvas, TString pdfFileName, bool &firstPage)
{
    // On the first page, ROOT requires adding "(" to open the file
    if (firstPage)
    {
        canvas->SaveAs(pdfFileName + "(");
        firstPage = false;
    }
    else
    {
        // For subsequent pages, save normally
        canvas->SaveAs(pdfFileName);
    }
}

///////////////////////////////////////////////////////////
// Function that visually marks peaks on the histogram canvas
///////////////////////////////////////////////////////////
void PaintPeaks(const std::vector<PeakInfo> &peaks_info, TH1 *hist)
{
    // Define vertical offsets for Y label placement (to avoid overlapping text)
    double y_offsets[] = {1.2, 1.1, 1.3, 1.05, 1.25, 1.15};
    int num_offsets = sizeof(y_offsets) / sizeof(y_offsets[0]);

    // Loop through each detected peak and annotate it on the histogram
    for (int i = 0; i < peaks_info.size(); ++i)
    {
        double x_peak = peaks_info[i].position;
        double y_peak_amplitude = peaks_info[i].amplitude;
        double y_label_pos = y_peak_amplitude * y_offsets[i % num_offsets]; // alternate label heights

        // Create a small text label showing the index of the peak
        TText *label = new TText(x_peak, y_label_pos, Form("%d", i));
        label->SetTextColor(kRed);
        label->SetTextFont(42);
        label->SetTextSize(0.03);
        label->Draw("same"); // draw on top of histogram
    }
}

///////////////////////////////////////////////////////////
// Function to detect the most relevant peaks in a histogram
// using ROOT's TSpectrum peak search
///////////////////////////////////////////////////////////
std::vector<PeakInfo> get_relevant_peaks(TSpectrum *spectrum, TH1 *hist,
                                         double sigma, double ratio_of_highest_peak,
                                         double c1, double c2)
{
    // Detect peaks automatically using TSpectrum
    Int_t nfound = spectrum->Search(hist, sigma, "", ratio_of_highest_peak);
    std::cout << "Number of peaks found: " << nfound << std::endl;

    if (nfound == 0)
        return {};

    double *xpeaks = spectrum->GetPositionX();
    double *ypeaks = spectrum->GetPositionY();

    std::vector<PeakInfo> peaks;

    // Analyze each detected peak and compute selection metrics
    for (int i = 0; i < nfound; ++i)
    {
        double x_peak = xpeaks[i];
        double y_peak = ypeaks[i];

        // Initial FWHM estimate assuming a Gaussian with given sigma
        double estimated_fwhm = 2.355 * sigma;

        // Estimate background level around the peak
        int bin_center = hist->FindBin(x_peak);
        int nBins = hist->GetNbinsX();

        int left = std::max(1, bin_center - int(c1 * sigma));
        int right = std::min(nBins, bin_center + int(c1 * sigma));

        double background_sum = 0.0;
        int background_count = 0;

        // Compute local background as average content around (but excluding) the peak bin
        for (int b = left; b <= right; ++b)
        {
            if (b == bin_center) continue;
            background_sum += hist->GetBinContent(b);
            background_count++;
        }

        double background_level = (background_count > 0) ? background_sum / background_count : 0.0;

        // Compute signal above background
        double signal = y_peak - background_level;

        // Estimate signal-to-noise ratio (simple sqrt-based model)
        double snr = (background_level > 0) ? signal / sqrt(background_level) : signal;

        std::cout << "Peak @ " << x_peak
                  << " | Amplitude: " << y_peak
                  << " | Bkg: " << background_level
                  << " | Signal: " << signal
                  << " | SNR: " << snr << std::endl;

        // Selection criterion: require sufficient signal and good SNR
        if (signal > c2 && snr > 9)
        {
            peaks.push_back({x_peak, estimated_fwhm, y_peak});
        }
    }

    // Sort peaks by amplitude (highest first)
    std::sort(peaks.begin(), peaks.end(),
              [](const auto &a, const auto &b) { return a.amplitude > b.amplitude; });

    return peaks;
}

///////////////////////////////////////////////////////////
// Function that checks whether an entry already exists in the output text file
// to avoid duplicate writing
///////////////////////////////////////////////////////////
bool entryExists(const TString &fileName, const TString &Source, int Detector_ID, double theoretical_energy)
{
    std::ifstream infile(fileName);
    if (!infile.is_open()) {
        return false; // file doesn’t exist yet
    }

    std::string line;
    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        std::string s_source;
        int s_id;
        double s_mean, s_energy, s_resolution;

        // Parse each line of the output text file
        if (!(iss >> s_source >> s_id >> s_mean >> s_energy >> s_resolution)) { continue; }

        // If we find the same source, detector, and theoretical energy — skip writing
        if (s_source == Source && s_id == Detector_ID && s_energy == theoretical_energy) {
            infile.close();
            return true;
        }
    }
    infile.close();
    return false;
}

///////////////////////////////////////////////////////////
// Function that fits a single peak with a Gaussian + linear background model
// and stores the results (mean, sigma, resolution) in a text file
///////////////////////////////////////////////////////////
void FitAndSaveToFile(TH1F *h, TString Detector, Int_t Detector_ID,
                      TString Source, Double_t theoretical_energy, TCanvas *canvas,
                      std::ofstream &outFile, std::vector<std::pair<double, double>> &res_points,
                      double peak_pos, double sigma_init, TString outputFolder)
{
    // Define the fitting window around the peak (proportional to peak position)
    float_t constant = 10.0 * peak_pos/2000;
    double x_min = peak_pos - constant * sigma_init;
    double x_max = peak_pos + constant * sigma_init;

    // Draw histogram and configure axis ranges
    canvas->cd();
    h->GetXaxis()->SetRangeUser(x_min * 0.8, x_max * 1.2);
    h->GetXaxis()->SetTitle("Channels");
    h->GetYaxis()->SetTitle("counts/bin");
    h->Draw();

    // Define a combined Gaussian + linear background function
    TF1 *gaus_expo = new TF1("gaus_expo_fit", "[0]*exp(-0.5*((x-[1])/[2])^2) + [3] + [4]*x", x_min, x_max);

    // Initialize parameters: amplitude, mean, sigma
    gaus_expo->SetParameter(0, h->GetBinContent(h->FindBin(peak_pos)));
    gaus_expo->SetParameter(1, peak_pos);
    gaus_expo->SetParameter(2, sigma_init);

    // Perform the fit (quiet mode, store result)
    TFitResultPtr fit = h->Fit(gaus_expo, "SQ+R", "", x_min, x_max);

    // Skip invalid fits
    if (!fit->IsValid())
    {
        std::cout << "*** Invalid fit for Detector " << Detector_ID << " at peak " << peak_pos << "\n";
        delete gaus_expo;
        return;
    }

    // Extract fitted parameters
    double mean = fit->Parameter(1);
    double sigma_fit = fit->Parameter(2);
    double FWHM = 2.355 * sigma_fit;
    double resolution = (FWHM / mean) * 100.0;

    // Save resolution vs energy point
    res_points.emplace_back(theoretical_energy, resolution);

    // Print detailed results to console
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "Detector " << Detector_ID << " | Mean: " << mean << " | Sigma: " << sigma_fit
              << " | FWHM: " << FWHM << " | Resolution: " << resolution << "%\n";

    // Prepare output filename for this detector
    TString outputFileName = Form("%sEdge_Compton_Detector=%s.txt", outputFolder.Data(), Detector.Data());

    // Avoid writing duplicate entries
    if (!entryExists(outputFileName, Source, Detector_ID, theoretical_energy)) {
        outFile << Source << " " << Detector_ID << " " << mean << " " << theoretical_energy << " " << resolution << std::endl;
    } else {
        std::cout << "Duplicate entry found for Detector " << Detector_ID << " and Source " << Source
                  << " with energy " << theoretical_energy << " keV. Skipping write." << std::endl;
    }

    // Draw results on canvas (text box with fit info)
    TPaveText *legend = new TPaveText(0.12, 0.7, 0.4, 0.88, "NDC");
    legend->SetBorderSize(0);
    legend->SetFillColor(0);
    legend->SetTextSize(0.03);
    legend->AddText(Form("Mean = %.2f", mean));
    legend->AddText(Form("Sigma = %.2f", sigma_fit));
    legend->AddText(Form("FWHM = %.2f", FWHM));
    legend->AddText(Form("Resolution = %.2f %%", resolution));
    legend->Draw("same");

    delete gaus_expo;
}

///////////////////////////////////////////////////////////
// Main calibration function
// Reads histograms, detects peaks, fits them, and saves results
///////////////////////////////////////////////////////////
void calibrations(TString FName, TString Detector, Int_t TotalDetector,
                  TString Source, const std::vector<int> &detector_ids,
                  const std::vector<double> &energy, double ratio_of_highest_peak,
                  double sigma_guess, double RangeUserXmin, double RangeUserXmax,
                  double c1, double c2, TString outputFolder)
{
    // Open input ROOT file containing histograms
    TFile *input_file = TFile::Open(FName);
    if (!input_file)
    {
        std::cout << "Error: Cannot open input file.\n";
        return;
    }

    // Open output text file in append mode (to keep previous runs)
    TString outputFileName = Form("%sEdge_Compton_Detector=%s.txt", outputFolder.Data(), Detector.Data());
    std::ofstream outFile;
    outFile.open(outputFileName, std::ios::app);
    if (!outFile.is_open()) {
        std::cerr << "Error: Could not open output file " << outputFileName << " for writing.\n";
        input_file->Close();
        return;
    }

    // Create objects for peak search and visualization
    TSpectrum *spectrum = new TSpectrum();
    TCanvas *canvas = new TCanvas("canvas", "Histogram", 800, 600);
    canvas->SetLogy(); // logarithmic Y scale to better see peaks

    // Prepare PDF file for visualization of all detectors
    TString pdfFileName = outputFolder + "calibration_results_" + Source + "_" + Detector + ".pdf";
    bool firstPage = true;

    // Set dry-run mode to true if only visualization is needed (no fitting)
    bool dry_run = false;

    // Data structure to store all energy-resolution pairs per detector
    std::map<int, std::vector<std::pair<double, double>>> all_resolutions;

    // Vector to hold all histograms read from file
    std::vector<TH1F *> histo_channels(TotalDetector + 1);

    // Loop through selected detector IDs
    for (auto i : detector_ids)
    {
        TString name = Form("histograma_Amplitude_%d", i);
        histo_channels[i] = (TH1F *)input_file->Get(name);
        if (!histo_channels[i])
        {
            std::cout << "Skipping detector " << i << " (no histogram found)\n";
            continue;
        }

        std::cout << "\n>>> Detector " << i << std::endl;

        // Clear and prepare canvas for each detector
        canvas->cd();
        canvas->Clear();
        histo_channels[i]->GetXaxis()->SetRangeUser(RangeUserXmin, RangeUserXmax);
        histo_channels[i]->Draw();

        // Detect candidate peaks using TSpectrum
        std::vector<PeakInfo> peaks_info = get_relevant_peaks(spectrum, histo_channels[i],
                                                              sigma_guess, ratio_of_highest_peak, c1, c2);

        // Skip detectors with no peaks
        if (peaks_info.empty())
        {
            std::cout << "No peaks found for detector " << i << ". Skipping.\n";
            SaveCanvasToPDF(canvas, pdfFileName, firstPage);
            continue;
        }

        // Annotate detected peaks on histogram
        PaintPeaks(peaks_info, histo_channels[i]);
        SaveCanvasToPDF(canvas, pdfFileName, firstPage);

        // If dry-run enabled, stop here (visual check only)
        if (dry_run)
        {
            std::cout << "Dry run complete. Please check the generated PDF file, and set `dry_run = false` "
                         "and `desired_peaks_per_detector` map with the correct indices before running again." << std::endl;
            continue;
        }

        // Prepare resolution vector for this detector
        all_resolutions[i].clear();

        // Fit each selected peak corresponding to known theoretical energies
        for (int p_index = 0; p_index < energy.size(); ++p_index)
        {
            if (p_index >= 0 && p_index < peaks_info.size())
            {
                // Special handling for AmBe source: use specific peak index
                if (Source == "AmBe") p_index = 2;

                const auto &peak = peaks_info[p_index];
                double peak_pos = peak.position;
                double sigma_for_fit = peak.fwhm / 2.355;
                double theoretical_energy = (p_index < energy.size()) ? energy[p_index] : 0.0;

                std::cout << "Fitting selected peak with index " << p_index << " at position " << peak_pos << std::endl;

                // Fit and save peak data to file
                FitAndSaveToFile(histo_channels[i], Detector, i, Source, theoretical_energy,
                                 canvas, outFile, all_resolutions[i], peak_pos, sigma_for_fit, outputFolder);

                // Save the fitted peak plot in the PDF
                SaveCanvasToPDF(canvas, pdfFileName, firstPage);
            }
            else
            {
                std::cout << "Invalid peak index: " << p_index << " for detector " << i << ". Skipping." << std::endl;
            }
        }
    }

    // Close PDF file properly
    canvas->Clear();
    canvas->SaveAs(pdfFileName + ")");

    // Close all open files and release memory
    outFile.close();
    delete canvas;
    delete spectrum;
    input_file->Close();

    cout << "\nCalibration data saved to " << outputFileName << std::endl;
}
