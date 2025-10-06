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

// O structură simplă pentru a stoca datele unui peak
struct PeakInfo
{
    double position;
    double fwhm;
    double amplitude;
};

// Funcție ajutătoare pentru a salva pe o pagină nouă în PDF
void SaveCanvasToPDF(TCanvas *canvas, TString pdfFileName, bool &firstPage)
{
    if (firstPage)
    {
        canvas->SaveAs(pdfFileName + "(");
        firstPage = false;
    }
    else
    {
        canvas->SaveAs(pdfFileName);
    }
}

// Funcția care desenează peak-urile pe canvas
void PaintPeaks(const std::vector<PeakInfo> &peaks_info, TH1 *hist)
{
    // Definește un set de offset-uri pentru poziția Y a etichetelor, pentru a evita suprapunerea
    double y_offsets[] = {1.2, 1.1, 1.3, 1.05, 1.25, 1.15};
    int num_offsets = sizeof(y_offsets) / sizeof(y_offsets[0]);

    for (int i = 0; i < peaks_info.size(); ++i)
    {
        double x_peak = peaks_info[i].position;
        double y_peak_amplitude = peaks_info[i].amplitude;
        double y_label_pos = y_peak_amplitude * y_offsets[i % num_offsets];

        TText *label = new TText(x_peak, y_label_pos, Form("%d", i));
        label->SetTextColor(kRed);
        label->SetTextFont(42);
        label->SetTextSize(0.03);
        label->Draw("same");
    }
}

// Funcție pentru a selecta cele mai proeminente peak-uri
std::vector<PeakInfo> get_relevant_peaks(TSpectrum *spectrum, TH1 *hist,
                                         double sigma, double ratio_of_highest_peak,
                                         double c1, double c2)
{
    Int_t nfound = spectrum->Search(hist, sigma, "", ratio_of_highest_peak);
    std::cout << "Number of peaks found: " << nfound << std::endl;

    if (nfound == 0)
        return {};

    double *xpeaks = spectrum->GetPositionX();
    double *ypeaks = spectrum->GetPositionY();

    std::vector<PeakInfo> peaks;

    for (int i = 0; i < nfound; ++i)
    {
        double x_peak = xpeaks[i];
        double y_peak = ypeaks[i];

        // Estimare FWHM inițial
        double estimated_fwhm = 2.355 * sigma;

        // Calculez background-ul ca media în jurul peak-ului
        int bin_center = hist->FindBin(x_peak);
        int nBins = hist->GetNbinsX();

        int left = std::max(1, bin_center - int(c1 * sigma));
        int right = std::min(nBins, bin_center + int(c1 * sigma));

        double background_sum = 0.0;
        int background_count = 0;
        for (int b = left; b <= right; ++b)
        {
            if (b == bin_center) continue; // să nu iau chiar binul cu peak-ul
            background_sum += hist->GetBinContent(b);
            background_count++;
        }
        double background_level = (background_count > 0) ? background_sum / background_count : 0.0;

        // Semnal peste fond
        double signal = y_peak - background_level;

        // Semnal-zgomot (simplu)
        double snr = (background_level > 0) ? signal / sqrt(background_level) : signal;

        std::cout << "Peak @ " << x_peak
                  << " | Amplitude: " << y_peak
                  << " | Bkg: " << background_level
                  << " | Signal: " << signal
                  << " | SNR: " << snr << std::endl;

        // Criteriu de selecție: semnal > c2 și SNR decent (>3 sigma)
        if (signal > c2 && snr > 9)
        {
            peaks.push_back({x_peak, estimated_fwhm, y_peak});
        }
    }

    std::sort(peaks.begin(), peaks.end(),
              [](const auto &a, const auto &b) { return a.amplitude > b.amplitude; });

    return peaks;
}


// === Funcție nouă pentru a verifica dacă o înregistrare există deja ===
bool entryExists(const TString &fileName, const TString &Source, int Detector_ID, double theoretical_energy)
{
    std::ifstream infile(fileName);
    if (!infile.is_open()) {
        return false;
    }

    std::string line;
    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        std::string s_source;
        int s_id;
        double s_mean, s_energy, s_resolution;

        if (!(iss >> s_source >> s_id >> s_mean >> s_energy >> s_resolution)) { continue; }

        if (s_source == Source && s_id == Detector_ID && s_energy == theoretical_energy) {
            infile.close();
            return true;
        }
    }
    infile.close();
    return false;
}

// Funcția de fitare a unui singur peak cu background exponențial, care salvează și datele
void FitAndSaveToFile(TH1F *h, TString Detector, Int_t Detector_ID,
                      TString Source, Double_t theoretical_energy, TCanvas *canvas,
                      std::ofstream &outFile, std::vector<std::pair<double, double>> &res_points, double peak_pos, double sigma_init, TString outputFolder)
{
    float_t constant = 10.0 * peak_pos/2000;
    double x_min = peak_pos - constant * sigma_init;
    double x_max = peak_pos + constant * sigma_init;

    canvas->cd();
    h->GetXaxis()->SetRangeUser(x_min * 0.8, x_max * 1.2);
    h->GetXaxis()->SetTitle("Channels");
    h->GetYaxis()->SetTitle("counts/bin");
    h->Draw();

    TF1 *gaus_expo = new TF1("gaus_expo_fit", "[0]*exp(-0.5*((x-[1])/[2])^2) + [3] + [4]*x", x_min, x_max);

    gaus_expo->SetParameter(0, h->GetBinContent(h->FindBin(peak_pos)));
    gaus_expo->SetParameter(1, peak_pos);
    gaus_expo->SetParameter(2, sigma_init);

    TFitResultPtr fit = h->Fit(gaus_expo, "SQ+R", "", x_min, x_max);

    if (!fit->IsValid())
    {
        std::cout << "*** Invalid fit for Detector " << Detector_ID << " at peak " << peak_pos << "\n";
        delete gaus_expo;
        return;
    }

    double mean = fit->Parameter(1);
    double sigma_fit = fit->Parameter(2);
    double FWHM = 2.355 * sigma_fit;
    double resolution = (FWHM / mean) * 100.0;

    res_points.emplace_back(theoretical_energy, resolution);

    std::cout << std::fixed << std::setprecision(2);
    std::cout << "Detector " << Detector_ID << " | Mean: " << mean << " | Sigma: " << sigma_fit
              << " | FWHM: " << FWHM << " | Resolution: " << resolution << "%\n";

    TString outputFileName = Form("%sEdge_Compton_Detector=%s.txt", outputFolder.Data(), Detector.Data());
    if (!entryExists(outputFileName, Source, Detector_ID, theoretical_energy)) {
        outFile << Source << " " << Detector_ID << " " << mean << " " << theoretical_energy << " " << resolution << std::endl;
    } else {
        std::cout << "Duplicate entry found for Detector " << Detector_ID << " and Source " << Source << " with energy " << theoretical_energy << " keV. Skipping write." << std::endl;
    }


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

void calibrations(TString FName, TString Detector, Int_t TotalDetector, TString Source, const std::vector<int> &detector_ids, const std::vector<double> &energy, double ratio_of_highest_peak, double sigma_guess, double RangeUserXmin, double RangeUserXmax, double c1, double c2, TString outputFolder)
{

    TFile *input_file = TFile::Open(FName);
    if (!input_file)
    {
        std::cout << "Error: Cannot open input file.\n";
        return;
    }

    // === Deschid fișierul o singură dată în modul append ===
    TString outputFileName = Form("%sEdge_Compton_Detector=%s.txt", outputFolder.Data(), Detector.Data());
    std::ofstream outFile;
    outFile.open(outputFileName, std::ios::app);
    if (!outFile.is_open()) {
        std::cerr << "Error: Could not open output file " << outputFileName << " for writing.\n";
        input_file->Close();
        return;
    }

    TSpectrum *spectrum = new TSpectrum();
    TCanvas *canvas = new TCanvas("canvas", "Histogram", 800, 600);
    canvas->SetLogy();

    TString pdfFileName = outputFolder + "calibration_results_" + Source + "_" + Detector + ".pdf";
    bool firstPage = true;

    bool dry_run = false;

    std::map<int, std::vector<std::pair<double, double>>> all_resolutions;
    std::vector<TH1F *> histo_channels(TotalDetector + 1);

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

        canvas->cd();
        canvas->Clear();
        histo_channels[i]->GetXaxis()->SetRangeUser(RangeUserXmin, RangeUserXmax);
        histo_channels[i]->Draw();

        std::vector<PeakInfo> peaks_info = get_relevant_peaks(spectrum, histo_channels[i], sigma_guess, ratio_of_highest_peak, c1, c2);

        if (peaks_info.empty())
        {
            std::cout << "No peaks found for detector " << i << ". Skipping.\n";
            SaveCanvasToPDF(canvas, pdfFileName, firstPage);
            continue;
        }

        PaintPeaks(peaks_info, histo_channels[i]);
        SaveCanvasToPDF(canvas, pdfFileName, firstPage);

        if (dry_run)
        {
            std::cout << "Dry run complete. Please check the generated PDF file, and set `dry_run = false` and `desired_peaks_per_detector` map with the correct indices before running again." << std::endl;
            continue;
        }

        all_resolutions[i].clear();

        for (int p_index = 0; p_index < energy.size(); ++p_index)
        {
            if (p_index >= 0 && p_index < peaks_info.size())
            {
                if (Source == "AmBe") p_index = 2;
                const auto &peak = peaks_info[p_index];
                double peak_pos = peak.position;
                double sigma_for_fit = peak.fwhm / 2.355;
                double theoretical_energy = (p_index < energy.size()) ? energy[p_index] : 0.0;

                std::cout << "Fitting selected peak with index " << p_index << " at position " << peak_pos << std::endl;
                FitAndSaveToFile(histo_channels[i], Detector, i, Source, theoretical_energy, canvas, outFile, all_resolutions[i], peak_pos, sigma_for_fit, outputFolder);

                SaveCanvasToPDF(canvas, pdfFileName, firstPage);
            }
            else
            {
                std::cout << "Invalid peak index: " << p_index << " for detector " << i << ". Skipping." << std::endl;
            }
        }
    }

    canvas->Clear();
    canvas->SaveAs(pdfFileName + ")");

    // === Închid fișierul la finalul scriptului ===
    outFile.close();

    delete canvas;
    delete spectrum;
    input_file->Close();

    cout << "\nCalibration data saved to " << outputFileName << std::endl;
}