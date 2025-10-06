#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TFitResult.h>
#include <TPaveText.h>
#include <TSpectrum.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TGraph.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <sstream>
#include <cmath> // Pentru TMath::Log și TMath::Power

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

// Structurile ajutătoare (nemodificate)
struct PeakInfo
{
    double position;
    double fwhm;
    double amplitude;
};

// !!! NOU: Structură pentru a ține datele sursei preluate din scriptul Bash/JSON
struct SourcePeakData
{
    std::string name;
    std::vector<double> energies;
    std::vector<std::pair<double, double>> search_windows;
    std::vector<double> sigma_guess;
    std::vector<double> intensity, dIntensity;
    double halflife, dHalflife, activity, dActivity;
    bool initial_activity = true;
    double measure_time, decay_time;
};

// Caută peak în fereastra dată cu TSpectrum (nemodificată)
PeakInfo FindPeakInWindow(TSpectrum *spectrum, TH1 *h, double sigma, double xmin, double xmax)
{
    h->GetXaxis()->SetRangeUser(xmin, xmax);

    double ratio_of_highest_peak = 0.05;
    Int_t nfound = spectrum->Search(h, sigma, "", ratio_of_highest_peak);
    if (nfound <= 0)
    {
        std::cout << "No peak found in range " << xmin << " - " << xmax << " keV\n";
        return {0, 0, 0};
    }

    double *xpeaks = spectrum->GetPositionX();
    double *ypeaks = spectrum->GetPositionY();

    // luam cel mai mare peak din fereastra
    int best = 0;
    for (int i = 1; i < nfound; i++)
    {
        if (ypeaks[i] > ypeaks[best])
            best = i;
    }

    double fwhm = 2.355 * sigma;
    return {xpeaks[best], fwhm, ypeaks[best]};
}

// Salvează canvas-ul într-un fișier PDF (nemodificată)
void SaveCanvasToPDF(TCanvas *canvas, const std::string &pdfFileName, bool &firstPage)
{
    if (firstPage)
    {
        canvas->Print((pdfFileName + "(").c_str());
        firstPage = false;
    }
    else
    {
        canvas->Print(pdfFileName.c_str());
    }
}

// Metoda alternativă de integrare a peak-ului cu background polinomial grad 1 (nemodificată)
void IntegratePeakWithBackground(TH1 *hist, double xmin_peak, double xmax_peak, double xmin_original, double xmax_original, TCanvas *canvas, std::ofstream &outFile)
{
    canvas->cd();
    hist->Draw();

    // --- Define background intervals ---
    int bin_xmin_bg_left = hist->FindBin(xmin_original);
    int bin_xmax_bg_left = hist->FindBin(xmin_peak);
    int bin_xmin_bg_right = hist->FindBin(xmax_peak);
    int bin_xmax_bg_right = hist->FindBin(xmax_original);

    // --- Collect points for background fit ---
    std::vector<double> x_bg, y_bg;
    for (int b = bin_xmin_bg_left; b <= bin_xmax_bg_left; ++b)
    {
        x_bg.push_back(hist->GetBinCenter(b));
        y_bg.push_back(hist->GetBinContent(b));
    }
    for (int b = bin_xmin_bg_right; b <= bin_xmax_bg_right; ++b)
    {
        x_bg.push_back(hist->GetBinCenter(b));
        y_bg.push_back(hist->GetBinContent(b));
    }

    // --- Fit polinom grad 1 pe background ---
    int npoints = x_bg.size();
    TGraph *gr_bg = new TGraph(npoints, &x_bg[0], &y_bg[0]);
    TF1 *bg_fit = new TF1("bg_fit", "[0]+[1]*x", xmin_original, xmax_original);
    TFitResultPtr r = gr_bg->Fit(bg_fit, "SQ"); // store fit result

    // --- Integrala sub peak ---
    int bin_xmin_peak = hist->FindBin(xmin_peak);
    int bin_xmax_peak = hist->FindBin(xmax_peak);

    double total_counts = 0;
    double bg_counts = 0;
    for (int b = bin_xmin_peak; b <= bin_xmax_peak; ++b)
    {
        double x = hist->GetBinCenter(b);
        double y = hist->GetBinContent(b);
        total_counts += y;

        double y_bg_est = bg_fit->Eval(x);
        bg_counts += y_bg_est;
    }

    double net_counts = total_counts - bg_counts;
    
    double total_counts_err = sqrt(total_counts);
    // double bg_counts_err    = sqrt(bg_counts);
    
    double dx  = xmax_peak - xmin_peak;
    double dx2 = 0.5*(xmax_peak*xmax_peak - xmin_peak*xmin_peak);
    double var_p0   = r->CovMatrix(0,0);
    double var_p1   = r->CovMatrix(1,1);
    double cov_p0p1 = r->CovMatrix(0,1);

    double bg_counts_err = sqrt( dx*dx*var_p0 + dx2*dx2*var_p1 + 2*dx*dx2*cov_p0p1 );

    double net_counts_err   = sqrt(total_counts_err*total_counts_err + bg_counts_err*bg_counts_err);


    // Scriere în fișier
    outFile << xmin_peak << ","
        << xmax_peak << ","
        << total_counts << ","
        << bg_counts << ","
        << net_counts << ","
        << total_counts_err << ","
        << bg_counts_err << ","
        << net_counts_err << ",";


    // --- Draw lines for visualization ---
    TLine *l1 = new TLine(xmin_peak, 0, xmin_peak, hist->GetMaximum());
    TLine *l2 = new TLine(xmax_peak, 0, xmax_peak, hist->GetMaximum());
    l1->SetLineColor(kRed);
    l1->SetLineStyle(2);
    l1->Draw();
    l2->SetLineColor(kRed);
    l2->SetLineStyle(2);
    l2->Draw();

    bg_fit->SetLineColor(kBlue);
    bg_fit->Draw("same");
}

// Returnează aria și eroarea acesteia (nemodificată)
std::vector<double> FitAndSaveToFile(TH1 *hist, const std::string &Detector, int id, const std::string &Source,
                                double theoretical_energy, TCanvas *canvas, std::ofstream &outFile,
                                std::vector<double> &resolutions, double mean_guess, double sigma_guess,
                                double xmin, double xmax, PeakInfo peak)
{
    canvas->cd();
    hist->Draw();

    // --- Gaussian + linear background ---
    TF1 *gaus_poly = new TF1("gaus_poly",
                             "[0]*exp(-0.5*((x-[1])/[2])^2) + [3] + [4]*x",
                             xmin, xmax);

    // Set initial parameters
    gaus_poly->SetParameter(0, hist->GetBinContent(hist->FindBin(mean_guess)));               // A
    gaus_poly->SetParameter(1, mean_guess);                                                   // mu
    gaus_poly->SetParameter(2, sigma_guess);                                                  // sigma
    gaus_poly->SetParLimits(2, 0, 1e6);                                                       // sigma >= 0
    gaus_poly->SetParameter(3, hist->GetBinContent(hist->FindBin(mean_guess - sigma_guess))); // p3 (const bg)
    gaus_poly->SetParameter(4, 0);                                                            // p4 (linear bg)

    // --- Fit ---
    TFitResultPtr fitResult = hist->Fit(gaus_poly, "RSQ"); // R=return, S=store, Q=quiet

    // --- Extract parameters ---
    double A     = gaus_poly->GetParameter(0);
    double mean  = gaus_poly->GetParameter(1);
    double sigma = gaus_poly->GetParameter(2);
    double p3    = gaus_poly->GetParameter(3);
    double p4    = gaus_poly->GetParameter(4);

    double FWHM = 2.355 * sigma;
    double resolution = (FWHM / mean) * 100.0;
    resolutions.push_back(resolution);

    // --- Integration window ---
    int constant = (FWHM < 45) ? 4 : 3;
    double integral_min = mean - constant*sigma;
    double integral_max = mean + constant*sigma;

    double binWidth = hist->GetXaxis()->GetBinWidth(1);

    // --- Integrals ---
    double total_integral = gaus_poly->Integral(integral_min, integral_max) / binWidth;

    // Background integral
    double bg_integral = ( p3*(integral_max - integral_min)
                       + 0.5*p4*(integral_max*integral_max - integral_min*integral_min) ) / binWidth;

    // Net area under Gaussian peak
    double area = total_integral - bg_integral;

    // --- Error propagation using covariance matrix ---
    // Partial derivatives
    double dA_area = sigma * sqrt(2*M_PI);
    double dsigma_area = A * sqrt(2*M_PI);

    // Covariance for Gaussian parameters
    double varA     = fitResult->CovMatrix(0,0);
    double varSigma = fitResult->CovMatrix(2,2);
    double covASigma= fitResult->CovMatrix(0,2);

    // Gaussian integral error
    double area_err_gauss = sqrt(dA_area*dA_area*varA
                                + dsigma_area*dsigma_area*varSigma
                                + 2*dA_area*dsigma_area*covASigma);

    // Background integral error
    double dx  = integral_max - integral_min;
    double dx2 = 0.5*(integral_max*integral_max - integral_min*integral_min);
    double varP3     = fitResult->CovMatrix(3,3);
    double varP4     = fitResult->CovMatrix(4,4);
    double covP3P4   = fitResult->CovMatrix(3,4);
    double area_err_bg = sqrt(dx*dx*varP3 + dx2*dx2*varP4 + 2*dx*dx2*covP3P4);

    // Total error on net area
    double area_err = sqrt(area_err_gauss*area_err_gauss + area_err_bg*area_err_bg);

    // --- Write to file ---
    outFile << Detector << ","
        << id << ","
        << Source << ","
        << theoretical_energy << ","
        << mean << ","
        << sigma << ","
        << FWHM << ","
        << resolution << ","
        << area << ","
        << area_err << ",";


    // --- Optional: parallel calculation ---
    IntegratePeakWithBackground(hist, integral_min, integral_max, xmin, xmax, canvas, outFile);

    // --- Draw TPaveText legend ---
    TPaveText *pt = new TPaveText(0.15, 0.75, 0.35, 0.88, "NDC");
    pt->SetTextSize(0.03);
    pt->SetFillStyle(0);
    pt->SetBorderSize(0);
    pt->AddText(Form("Mean = %.2f", mean));
    pt->AddText(Form("FWHM = %.2f", FWHM));
    pt->AddText(Form("Res = %.2f %%", resolution));
    pt->AddText(Form("Area = %.1f +/- %.1f", area, area_err));
    pt->Draw();

    return {area, area_err};
}


// Calculează eficiența cu propagarea erorilor (nemodificată)
void efficiency_calc(double halflife, double dHalflife, double activity, double dActivity, 
                     double intensity, double dIntensity, 
                     double area, double dArea, std::ofstream &outFile, 
                     bool initial_activity, double measure_time, double decay_time)
{
    double lambda = log(2) / halflife;

    // activitatea la momentul masurarii
    double activity_exp;
    if(initial_activity)
        activity_exp = activity * exp(-lambda * decay_time); // activitatea finala
    else
        activity_exp = activity * exp(lambda * decay_time); // activitatea initiala

    double I_gamma = intensity / 100.0;
    double efficiency = area / (activity_exp * I_gamma * measure_time);

    // -------------------------------
    // Propagarea erorilor
    // -------------------------------
    double rel_err_area = dArea / area;

    double rel_err_activity = (dActivity > 0) ? (dActivity / activity) : 0.0;
    double dLambda = (log(2) / (halflife * halflife)) * dHalflife;
    double rel_err_T12 = decay_time * dLambda;
    double rel_err_activity_exp = sqrt(rel_err_activity*rel_err_activity + rel_err_T12*rel_err_T12);
    double rel_err_Igamma   = (dIntensity / intensity);

    // eroare relativa totala
    double rel_err_total = sqrt(pow(rel_err_area, 2) + pow(rel_err_activity_exp, 2) + pow(rel_err_Igamma, 2));

    double eff_err = efficiency * rel_err_total;

    // -------------------------------
    // Afisare
    outFile << halflife << ","
        << activity << ","
        << dActivity << ","
        << activity_exp << ","
        << activity_exp * rel_err_activity_exp << ","
        << intensity << ","
        << dIntensity << ","
        << measure_time << ","
        << efficiency * 100 << ","
        << eff_err * 100 << std::endl;

}

// Analiza peak-urilor dintr-un fișier ROOT
// !!! NOU: Am adăugat argumentele pentru SourcePeakData
void peak_analysis(const std::string &rootFileName, const std::string &Detector, const std::string &Source, const std::vector<int> &ids, const SourcePeakData &sourceData)
{
    TFile *input_file = TFile::Open(rootFileName.c_str(), "READ");
    if (!input_file || input_file->IsZombie())
    {
        std::cerr << "Error opening file " << rootFileName << std::endl;
        return;
    }

    std::string folder = "efficiency/";
    createDirectories(folder.c_str());
    //mkdir(folder.c_str(), 0777);

    std::map<int, TH1F *> histo_channels;
    std::map<int, std::vector<double>> all_resolutions;
    TSpectrum *spectrum = new TSpectrum();

    // !!! VECHE: Era hardcodarea datelor surselor aici. NOU: Le primim direct prin 'sourceData'.
    
    // Verificare rapidă a datelor
    if (sourceData.energies.empty() || sourceData.energies.size() != sourceData.search_windows.size()) {
        std::cerr << "Eroare: Datele sursei sunt invalide sau lipsesc (Energies/Search Windows)." << std::endl;
        input_file->Close();
        return;
    }

    std::string outputFileName = folder + "PeakAreas_Detector=" + Detector + "_Source=" + Source + ".txt";
    std::ofstream outFile(outputFileName);
    std::string pdfFileName = folder + "Fits_" + Source + "_" + Detector + ".pdf";
    TCanvas *canvas = new TCanvas("canvas", "Fits", 800, 600);
    bool firstPage = true;

    // Antetul fișierului de ieșire
    outFile << "Detector,ID,Source,Energy,Mean,Sigma,FWHM,Resolution,Area,AreaErr,"
                    << "PeakIntervalMin,PeakIntervalMax,TotalCounts,BgCounts,NetCounts,"
                    << "TotalCountsErr,BgCountsErr,NetCountsErr,"
                    << "Halflife(s),InitialActivity(Bq),dActivity,ActivityAtMeas(Bq),dActivityAtMeas,"
                    << "Intensity(%),dIntensity,MeasureTime(s),Efficiency(%),EfficiencyErr(%)"
                    << std::endl;

    for (auto i : ids)
    {
        // Presupunem că histograma este deja calibrată în energie (de la Calibrarea 4)
        TString name = Form("histograma_Energy_%d", i);
        histo_channels[i] = (TH1F *)input_file->Get(name);
        if (!histo_channels[i])
        {
            std::cout << "Missing histogram for det " << i << ". Continuăm." << std::endl;
            continue;
        }

        std::cout << "\n>>> Detector " << i << std::endl;
        all_resolutions[i].clear();

        for (int p = 0; p < sourceData.search_windows.size(); ++p)
        {
            auto [xmin, xmax] = sourceData.search_windows[p];
            // Folosim peak-ul găsit ca estimare de bază
            PeakInfo peak = FindPeakInWindow(spectrum, histo_channels[i], sourceData.sigma_guess[p], xmin, xmax);
            if (peak.amplitude <= 0)
            {
                std::cout << "Nu s-a găsit peak in intervalul " << xmin << " - " << xmax << " keV pentru det " << i << ". Sarim." << std::endl;
                continue;
            }

            double sigma_for_fit = peak.fwhm / 2.355;
            double theoretical_energy = sourceData.energies[p];

            // 1. Fit și calcul al ariei
            std::vector<double> area = FitAndSaveToFile(histo_channels[i], Detector, i, Source,
                             theoretical_energy, canvas, outFile,
                             all_resolutions[i], peak.position, sigma_for_fit, xmin, xmax, peak);

            // 2. Salvarea graficului
            SaveCanvasToPDF(canvas, pdfFileName, firstPage);

            // 3. Calculul eficienței
            efficiency_calc(sourceData.halflife, sourceData.dHalflife, sourceData.activity, sourceData.dActivity, sourceData.intensity[p],
                sourceData.dIntensity[p], area[0], area[1], outFile, sourceData.initial_activity, sourceData.measure_time, sourceData.decay_time);
        }
    }

    canvas->Print((pdfFileName + "]").c_str()); // închidem PDF
    input_file->Close();

    std::cout << "✅ Rezultatele au fost salvate in: " << outputFileName << std::endl;
}