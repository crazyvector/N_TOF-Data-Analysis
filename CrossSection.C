// ===========================================================
// ROOT Macro: Extract cross section for a gamma line of interest
// ===========================================================
#include <TH2D.h>
#include <TH1D.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TSpectrum.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <iostream>
#include <vector>
#include <sys/stat.h>
#include <sys/types.h>
#include <sstream>

// ===========================================================
// Helper function: Recursively create directories if they don't exist
// ===========================================================
void createDirectories(const std::string &path) {
    std::stringstream ss(path);
    std::string item;
    std::string currentPath = "."; // start from current directory

    while (std::getline(ss, item, '/')) {
        if (item.empty()) continue;
        currentPath += "/" + item;
        mkdir(currentPath.c_str(), 0777); // ignore error if folder exists
    }
}

// ===========================================================
// Undo log binning for flux histogram
// Converts histogram content to physical units using number of dedicated pulses
// ===========================================================
void UndoBinning(TH1D *hin, Double_t nr_dedicated_pulses)
{
    for (int j = 1; j <= hin->GetNbinsX(); j++)
    {
        Double_t bc = hin->GetBinContent(j);
        Double_t be = hin->GetBinError(j);
        Double_t bl = hin->GetBinLowEdge(j);
        Double_t bu = bl + hin->GetBinWidth(j);
        Double_t logbw = TMath::Log(bu / bl);
        hin->SetBinContent(j, bc * logbw * nr_dedicated_pulses);
        hin->SetBinError(j, be * logbw * nr_dedicated_pulses);
    }
}

// ===========================================================
// Save a canvas to PDF, handling multiple pages
// ===========================================================
void SaveCanvasToPDF(TCanvas *canvas, TString pdfFileName, bool &firstPage)
{
    if (firstPage)
    {
        canvas->SaveAs(pdfFileName + "("); // first page
        firstPage = false;
    }
    else
    {
        canvas->SaveAs(pdfFileName); // append subsequent pages
    }
}

// ===========================================================
// Main function: Compute cross section for interest gamma line
// ===========================================================
void CrossSection(TString inputFileName, TString histogramName,
    TString inputFluxName, TString histogramFluxName,
    double interestGammaEnergy, int rebinFactor,
    int minEntries, double nr_dedicated_pulses,
    double massNr, double rho,
    double eff, double peakSigma,
    double peakRatio, double preFitRange,
    double gammaWindow, TString pdfName,
    TString outputRootName)
{
    // --- Open input files ---
    TFile *f = new TFile(inputFileName, "READ");
    TH2D *h2 = (TH2D *)f->Get(histogramName); // neutron energy vs gamma energy
    if (!h2) { 
        std::cerr << "Histogram " << histogramName << " not found in file!" << std::endl; 
        return; 
    }

    TFile *f2 = new TFile(inputFluxName, "READ");
    TH1D *hFlux = (TH1D *)f2->Get(histogramFluxName);
    if (!hFlux) { 
        std::cerr << "Histogram " << histogramFluxName << " not found in file!" << std::endl;
        return;
    }

    int nBinsX = h2->GetNbinsX();

    // --- Results containers ---
    std::vector<double> E_n, sigma, Eerr, sigmaErr;
    bool firstPage = true;

    // --- Loop over neutron energy bins ---
    for (int ix = 1; ix <= nBinsX; ix += rebinFactor)
    {
        // neutron energy in eV
        //gamma energy in KeV
        //interest gamma energy in KeV
        double E_n_low = h2->GetXaxis()->GetBinLowEdge(ix);
        double E_n_high = h2->GetXaxis()->GetBinUpEdge(std::min(ix + rebinFactor - 1, nBinsX));
        double Ecenter = 0.5 * (E_n_low + E_n_high);
        double bin_width = E_n_high - E_n_low;

        std::cout << "Neutron Energy = " << Ecenter / 1e3 << " KeV"
                  << ", interest gamma energy = " << interestGammaEnergy << " KeV" << std::endl;

        if (Ecenter < interestGammaEnergy * 1e3) continue;

        int superBinNum = (ix - 1) / rebinFactor + 1;
        int nSuperBins = (nBinsX - 1) / rebinFactor + 1;
        std::cout << "Processing super-bin " << superBinNum << "/" << nSuperBins
                  << " (bins " << ix << "-" << std::min(ix + rebinFactor - 1, nBinsX) << ")" << std::endl;

        // --- Project gamma spectrum for this super-bin ---
        TH1D *projY = h2->ProjectionY(Form("projY_%d", ix), ix, std::min(ix + rebinFactor - 1, nBinsX));
        projY->GetXaxis()->SetRangeUser(interestGammaEnergy * (1.0 - 0.5), interestGammaEnergy * (1.0 + 0.5));

        if (projY->GetEntries() < minEntries) { // skip empty bins
            delete projY;
            continue;
        }

        double peakE = interestGammaEnergy; // assume known peak for now

        // =====================================================
        // Two-step Gaussian fit with safety checks
        // =====================================================

        // Step 1: Preliminary fit
        double fitLow_pre = peakE - preFitRange;
        double fitHigh_pre = peakE + preFitRange;

        TF1 *fitPre = new TF1("fitPre",
                              "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]+[4]*x",
                              fitLow_pre, fitHigh_pre);
        fitPre->SetParameters(projY->GetMaximum(), peakE, 50, projY->GetMinimum(), 0.0);
        
        TFitResultPtr fitResPre = projY->Fit(fitPre, "SQ0R", "", fitLow_pre, fitHigh_pre);
        if (fitResPre.Get() == nullptr || !fitResPre->IsValid()) {
            std::cout << "Pre-fit failed. Skipping bin." << std::endl;
            delete projY;
            continue;
        }

        double mean_pre = fitResPre->Parameter(1);
        double sigma_pre = fitResPre->Parameter(2);

        double xMinHist = projY->GetXaxis()->GetXmin();
        double xMaxHist = projY->GetXaxis()->GetXmax();
        double fitLow = std::max(xMinHist, mean_pre - 3 * sigma_pre);
        double fitHigh = std::min(xMaxHist, mean_pre + 3 * sigma_pre);
        cout << "Fit range: [" << fitLow << ", " << fitHigh << "]" << endl;

        if (fitHigh <= fitLow) {
            std::cout << "Invalid fit range. Skipping bin." << std::endl;
            delete projY;
            continue;
        }

        // Step 2: Final fit
        TF1 *fitFunc = new TF1("fitFunc",
                               "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]+[4]*x",
                               fitLow, fitHigh);
        fitFunc->SetParameters(fitResPre->Parameter(0), mean_pre, sigma_pre,
                               fitResPre->Parameter(3), fitResPre->Parameter(4));

        TFitResultPtr fit = projY->Fit(fitFunc, "SQ+R", "", fitLow, fitHigh);
        if (fit.Get() == nullptr || !fit->IsValid()) {
            std::cout << "Final fit failed. Skipping bin." << std::endl;
            delete projY;
            continue;
        }

        // --- Extract Gaussian parameters ---
        double A = fit->Parameter(0);
        double mean = fit->Parameter(1);
        double sigmaG = fit->Parameter(2);

        // --- Gaussian area and uncertainty ---
        const double SQRT2PI = sqrt(2.0 * TMath::Pi());
        double binwidth = projY->GetBinWidth(1); // assume uniform bin width
        double area = A * sigmaG * SQRT2PI / binwidth;
        // double area = fitFunc->Integral(fitLow, fitHigh);

        double varA = pow(fit->ParError(0), 2);
        double varS = pow(fit->ParError(2), 2);
        double covAS = 0.0;
        try { covAS = fit->CovMatrix(0, 2); } catch (...) { covAS = 0.0; }

        double dA = sigmaG * SQRT2PI;
        double dS = A * SQRT2PI;
        double areaVar = dA*dA*varA + dS*dS*varS + 2.0*dA*dS*covAS;
        double areaErr = sqrt(std::max(areaVar, 0.0));

        std::cout << "Peak final: " << mean << " ± " << fit->ParError(1) << " keV" << std::endl;
        std::cout << "Gaussian area = " << area << " ± " << areaErr << std::endl;

        // --- Integrate neutron flux for this super-bin ---
        int binLow = hFlux->FindBin(E_n_low);
        int binHigh = hFlux->FindBin(E_n_high);

        // Integrate flux using bin widths (ROOT handles non-uniform bins correctly)
        double neutronFlux = hFlux->Integral(binLow, binHigh, "width"); // use GetBinContent

        // Compute variance of flux (error propagation)
        double fluxVar = 0.0;
        for (int jx = binLow; jx <= binHigh; ++jx) {
            double err = hFlux->GetBinError(jx);
            double width = hFlux->GetBinWidth(jx);
            fluxVar += pow(err * width, 2);
        }

        double neutronFluxErr = (neutronFlux > 0) ? sqrt(fluxVar) : 0.0;

        if (neutronFlux <= 0) {
            std::cout << "Neutron flux zero. Skipping bin." << std::endl;
            delete projY;
            continue;
        }

        // --- Cross section constants ---
        const double amu = 1.66e-24;  // atomic mass unit [g]
        const double barn = 1e-24;    // cm^2
        double CONST = (amu * massNr) / (eff * rho) * (1.0 / (4.0 * TMath::Pi()));  // cm^2
        CONST /= barn; // convert to barns

        // --- Cross section calculation ---
        double xs = (area / neutronFlux) * CONST * 1000;

        // Statistical uncertainty propagation
        double xsVar_stat = pow(CONST, 2) * (
            areaVar / pow(neutronFlux, 2) +
            pow(area, 2) * fluxVar / pow(neutronFlux, 4)
        );

        if (xsVar_stat < 0) xsVar_stat = 0;
        double xsErr_stat = sqrt(xsVar_stat);

        // Optional: systematic uncertainties (currently zero)
        double relEffErr = 0.0;
        double relRhoErr = 0.0;
        double xsErr_sys = xs * sqrt(relEffErr * relEffErr + relRhoErr * relRhoErr);

        // Total uncertainty
        double xsErr_total = sqrt(xsErr_stat * xsErr_stat + xsErr_sys * xsErr_sys);

        // --- Store results ---
        E_n.push_back(Ecenter);
        sigma.push_back(xs);
        Eerr.push_back(0.5*bin_width);
        sigmaErr.push_back(xsErr_total);

        projY->GetXaxis()->SetRangeUser(fitLow, fitHigh);

        // --- Display progress ---
        projY->Draw();
        fitFunc->Draw("same");

        TPaveText *pt = new TPaveText(0.15, 0.7, 0.5, 0.9, "NDC");
        pt->SetTextSize(0.03);
        pt->SetFillColor(0);
        pt->SetBorderSize(1);
        pt->AddText(Form("Mean = %.2f ± %.2f keV", mean, fit->ParError(1)));
        pt->AddText(Form("Sigma = %.2f ± %.2f keV", sigmaG, fit->ParError(2)));
        pt->AddText(Form("Area = %.1f ± %.1f", area, areaErr));
        pt->Draw();

        gPad->Update();
        SaveCanvasToPDF((TCanvas*)gPad, pdfName, firstPage);

        std::cout << "Area = " << area << ", Flux = " << neutronFlux
                  << ", xs = " << xs << " ± " << xsErr_total << " Const: " << CONST << std::endl;

        delete projY;
        std::cout << std::endl;
    }

    // --- Final graph: cross section vs neutron energy ---
    int N = E_n.size();
    TGraphErrors *gr = new TGraphErrors(N, &E_n[0], &sigma[0], &Eerr[0], &sigmaErr[0]);
    TString title = Form("Cross section for %.1f keV gamma; Neutron energy [eV]; Cross Section [mb]", interestGammaEnergy);
    gr->SetTitle(title);
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(kBlue);

    TCanvas *c1 = new TCanvas("c1", "Cross section", 800, 600);
    gr->Draw("AP");
    c1->SaveAs(pdfName + ")");

    // --- Save output ROOT file ---
    TFile *fout = new TFile(outputRootName, "RECREATE");
    gr->SetName("cross_section");
    gr->Write("sigma");
    fout->Close();

    std::cout << "Analysis done. Graph saved to " << outputRootName << std::endl;
}

// ===========================================================
// Helper function: Load and correct flux histogram
// ===========================================================
void RunGetFlux(TString inputFluxName, TString histogramName, TString outputFluxName, double nr_dedicated_pulses)
{
    try {
        TFile *f = new TFile(inputFluxName, "READ");
        TH1D *h = (TH1D *)f->Get(histogramName);
        if (!h) { std::cerr << "Histogram " << histogramName << " not found in file!" << std::endl; f->Close(); return; }

        TH1D *hNew = (TH1D *)h->Clone(histogramName.Data());
        UndoBinning(hNew, nr_dedicated_pulses);

        TFile *fout = new TFile(outputFluxName, "RECREATE");
        hNew->SetName(histogramName);
        hNew->Write();
        fout->Close();
        f->Close();

        std::cout << "Flux histogram written to " << outputFluxName << " successfully." << std::endl;
    }
    catch (const std::exception &e) {
        std::cerr << "An exception occurred: " << e.what() << std::endl;
    }
    catch (...) {
        std::cerr << "An unknown error occurred." << std::endl;
    }
}

// ===========================================================
// Wrapper function: run flux correction and cross section calculation
// ===========================================================
void RunCrossSection(TString inputFileName,
    TString histogramName,
    TString inputFluxName,
    TString outputFluxName,
    TString histogramFluxName,
    double interestGammaEnergy,
    int rebinFactor,
    int minEntries,
    double nr_dedicated_pulses,
    double massNr,
    double density,
    double eff,
    double peakSigma,
    double peakRatio,
    double preFitRange,
    double gammaWindow,
    TString pdfName,
    TString outputRootName)
{
    try {
        std::string folder = "cross_section/";
        createDirectories(folder.c_str());

        RunGetFlux(inputFluxName, histogramFluxName, folder + outputFluxName, nr_dedicated_pulses);

        CrossSection(inputFileName, histogramName,
                     inputFluxName, histogramFluxName,
                     interestGammaEnergy, rebinFactor,
                     minEntries, nr_dedicated_pulses,
                     massNr, density, eff, peakSigma, peakRatio, preFitRange, gammaWindow,
                     folder + pdfName, folder + outputRootName);
    }
    catch (const std::exception &e) {
        std::cerr << "An exception occurred: " << e.what() << std::endl;
    }
    catch (...) {
        std::cerr << "An unknown error occurred." << std::endl;
    }
}
