#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TFitResult.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TGraph.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <sstream>
#include <cmath>
#include <sys/stat.h>
#include <sys/types.h>

// Recursively creates directories based on a given path
void createDirectories(const std::string &path) {
    std::stringstream ss(path);
    std::string item;
    std::string currentPath = ".";
    while (std::getline(ss, item, '/')) {
        if (item.empty()) continue;
        currentPath += "/" + item;
        mkdir(currentPath.c_str(), 0777);
    }
}

// Structure holding source data (energy, windows, intensities, activity, etc.)
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

// Saves the canvas to a PDF file
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

// Integrates the area under a peak using a first-degree polynomial background
// (Extended comments added; code remains identical)
void IntegratePeakWithBackground(TH1 *hist, double xmin_peak, double xmax_peak, double xmin_original, double xmax_original, TCanvas *canvas, std::ofstream &outFile)
{
    // Draw the histogram on the canvas (not used for numeric calculation, only for visualization)
    canvas->cd();
    hist->Draw();

    int bin_xmin_bg_left = hist->FindBin(xmin_original);
    int bin_xmax_bg_left = hist->FindBin(xmin_peak);
    int bin_xmin_bg_right = hist->FindBin(xmax_peak);
    int bin_xmax_bg_right = hist->FindBin(xmax_original);

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

    // Check if there are non-zero background data
    bool has_bg_data = false;
    for (double y : y_bg) {
        if (y > 0) { has_bg_data = true; break; }
    }
    if (!has_bg_data) {
        std::cout << "Warning: Background data is empty in the interval " << xmin_original << " - " << xmax_original << ".\n";
        outFile << xmin_peak << "," << xmax_peak << ",NaN,NaN,NaN,NaN,NaN,NaN,";
        return;
    }

    int npoints = x_bg.size();
    TGraph *gr_bg = new TGraph(npoints, &x_bg[0], &y_bg[0]);
    TF1 *bg_fit = new TF1("bg_fit", "[0]+[1]*x", xmin_original, xmax_original);
    TFitResultPtr r = gr_bg->Fit(bg_fit, "SQ");

    int bin_xmin_peak = hist->FindBin(xmin_peak);
    int bin_xmax_peak = hist->FindBin(xmax_peak);

    double total_counts = 0, bg_counts = 0;
    for (int b = bin_xmin_peak; b <= bin_xmax_peak; ++b)
    {
        double x = hist->GetBinCenter(b);
        double y = hist->GetBinContent(b);
        total_counts += y;
        double y_bg_est = bg_fit->Eval(x);
        bg_counts += y_bg_est;
    }

    // net_counts = total_counts - bg_counts
    double net_counts = total_counts - bg_counts;

    // Error for total_counts: assume Poisson -> sigma = sqrt(N)
    // (note: for small counts or preprocessing, another estimate can be used)
    double total_counts_err = sqrt(total_counts);

    // To calculate the error of the background estimate integrated over the peak interval,
    // the expression for the area under the background line is noted:
    //   S_bg = p0 * dx + p1 * (1/2)*(x2^2 - x1^2)
    // For error propagation, we need derivatives with respect to p0 and p1.
    double dx  = xmax_peak - xmin_peak;
    double dx2 = 0.5*(xmax_peak*xmax_peak - xmin_peak*xmin_peak);

    // Check valid pointer for fitResult
    if (!r.Get()) {
        outFile << xmin_peak << "," << xmax_peak << "," << total_counts << "," << bg_counts << "," << net_counts << ","
                << total_counts_err << ",NaN,NaN,";
        return;
    }

    // Extract covariance matrix elements from the background fit
    double var_p0   = r->CovMatrix(0,0); // Var(p0)
    double var_p1   = r->CovMatrix(1,1); // Var(p1)
    double cov_p0p1 = r->CovMatrix(0,1); // Cov(p0,p1)

    // Error of background area: sqrt( (dx)^2 Var(p0) + (dx2)^2 Var(p1) + 2 dx dx2 Cov(p0,p1) )
    // Derivative S_bg / dp0 = dx,    S_bg / dp1 = dx2
    double bg_counts_err = sqrt( dx*dx*var_p0 + dx2*dx2*var_p1 + 2*dx*dx2*cov_p0p1 );

    // Error of net_counts combines numerator error (Poisson) with background error (assumed independent)
    double net_counts_err   = sqrt(total_counts_err*total_counts_err + bg_counts_err*bg_counts_err);

    outFile << xmin_peak << ","
        << xmax_peak << ","
        << total_counts << ","
        << bg_counts << ","
        << net_counts << ","
        << total_counts_err << ","
        << bg_counts_err << ","
        << net_counts_err << ",";

    TLine *l1 = new TLine(xmin_peak, 0, xmin_peak, hist->GetMaximum());
    TLine *l2 = new TLine(xmax_peak, 0, xmax_peak, hist->GetMaximum());
    l1->SetLineColor(kRed); l1->SetLineStyle(2); l1->Draw();
    l2->SetLineColor(kRed); l2->SetLineStyle(2); l2->Draw();
    bg_fit->SetLineColor(kBlue); bg_fit->Draw("same");
}

std::vector<double> FitAndSaveToFile(TH1 *hist, const std::string &Detector, int id, const std::string &Source,
                                double theoretical_energy, TCanvas *canvas, std::ofstream &outFile,
                                std::vector<double> &resolutions, double mean_guess, double sigma_guess,
                                double xmin, double xmax)
{
    canvas->cd();
    hist->GetXaxis()->SetRangeUser(xmin, xmax);
    hist->Draw();

    // Check if there are non-zero data in the fit interval
    int bin_min = hist->FindBin(xmin);
    int bin_max = hist->FindBin(xmax);
    bool has_data = false;
    for (int b = bin_min; b <= bin_max; ++b) {
        if (hist->GetBinContent(b) > 0) {
            has_data = true;
            break;
        }
    }
    if (!has_data) {
        std::cout << "Warning: Histogram is empty in the interval " << xmin << " - " << xmax << ". Fit will be skipped.\n";
        outFile << Detector << "," << id << "," << Source << "," << theoretical_energy << ","
                << "NaN,NaN,NaN,NaN,NaN,NaN,";
        return {0, 0};
    }

    // Define Gaussian function + first-degree polynomial for background:
    //   f(x) = A * exp(-0.5 * ((x - mean)/sigma)^2) + p3 + p4 * x
    // Parameters: [0]=A, [1]=mean, [2]=sigma, [3]=p3, [4]=p4
    TF1 *gaus_poly = new TF1("gaus_poly",
                             "[0]*exp(-0.5*((x-[1])/[2])^2) + [3] + [4]*x",
                             xmin, xmax);

    gaus_poly->SetParameter(0, hist->GetBinContent(hist->FindBin(mean_guess)));
    gaus_poly->SetParameter(1, mean_guess);
    gaus_poly->SetParameter(2, sigma_guess);
    gaus_poly->SetParLimits(2, 0, 1e6);
    gaus_poly->SetParameter(3, hist->GetBinContent(hist->FindBin(mean_guess - sigma_guess)));
    gaus_poly->SetParameter(4, 0);

    TFitResultPtr fitResult = hist->Fit(gaus_poly, "RSQ");

    // Check valid pointer for fitResult
    if (!fitResult.Get()) {
        std::cout << "Warning: Fit did not return valid result for detector " << id << " at energy " << theoretical_energy << ".\n";
        outFile << Detector << "," << id << "," << Source << "," << theoretical_energy << ","
                << "NaN,NaN,NaN,NaN,NaN,NaN,";
        return {0, 0};
    }

    // Extract fit parameters (same notation as original code)
    double A     = gaus_poly->GetParameter(0);
    double mean  = gaus_poly->GetParameter(1);
    double sigma = gaus_poly->GetParameter(2);
    double p3    = gaus_poly->GetParameter(3);
    double p4    = gaus_poly->GetParameter(4);

    double FWHM = 2.355 * sigma;
    double resolution = (FWHM / mean) * 100.0;
    resolutions.push_back(resolution);

    int constant = (FWHM < 45) ? 4 : 3;
    double integral_min = mean - constant*sigma;
    double integral_max = mean + constant*sigma;
    double binWidth = hist->GetXaxis()->GetBinWidth(1);

    // Numeric calculation of the area (integral of the continuous function divided by bin width)
    double total_integral = gaus_poly->Integral(integral_min, integral_max) / binWidth;

    // Area under background (linear polynomial):
    //   S_bg = p3*(x2-x1) + 0.5*p4*(x2^2 - x1^2)
    double bg_integral = ( p3*(integral_max - integral_min)
                       + 0.5*p4*(integral_max*integral_max - integral_min*integral_min) ) / binWidth;
    double area = total_integral - bg_integral;

    // ============================
    // GAUSSIAN AREA ERROR
    // ============================
    // Theoretical formula: for Gaussian f(x) = A * exp(-((x-mean)^2)/(2 sigma^2))
    // Full area: S = A * sigma * sqrt(2*pi)
    // Derivatives:
    //   dS/dA = sigma * sqrt(2*pi)
    //   dS/dsigma = A * sqrt(2*pi)
    // Propagation (including covariance between A and sigma):
    //   (ΔS)^2 = (dS/dA)^2 Var(A) + (dS/dsigma)^2 Var(sigma) + 2 (dS/dA)(dS/dsigma) Cov(A,sigma)
    // Implementation (use covariance matrix elements from fitResult):
    double dA_area = sigma * sqrt(2*M_PI);
    double dsigma_area = A * sqrt(2*M_PI);
    double varA     = fitResult->CovMatrix(0,0); // Var(A)
    double varSigma = fitResult->CovMatrix(2,2); // Var(sigma)
    double covASigma= fitResult->CovMatrix(0,2); // Cov(A,sigma)
    double area_err_gauss = sqrt(dA_area*dA_area*varA
                                + dsigma_area*dsigma_area*varSigma
                                + 2*dA_area*dsigma_area*covASigma);

    // ============================
    // BACKGROUND AREA ERROR (LINEAR POLYNOMIAL)
    // ============================
    // Reminder: S_bg = p3*dx + p4*(1/2)*(x2^2 - x1^2)
    // Derivatives:
    //   dS_bg/dp3 = dx
    //   dS_bg/dp4 = 0.5*(x2^2 - x1^2) = dx2
    // Propagation:
    //   (ΔS_bg)^2 = (dx)^2 Var(p3) + (dx2)^2 Var(p4) + 2 dx dx2 Cov(p3,p4)
    double dx  = integral_max - integral_min;
    double dx2 = 0.5*(integral_max*integral_max - integral_min*integral_min);
    double varP3     = fitResult->CovMatrix(3,3);
    double varP4     = fitResult->CovMatrix(4,4);
    double covP3P4   = fitResult->CovMatrix(3,4);
    double area_err_bg = sqrt(dx*dx*varP3 + dx2*dx2*varP4 + 2*dx*dx2*covP3P4);

    // Total area error (gauss + background) combined as sum of variances
    double area_err = sqrt(area_err_gauss*area_err_gauss + area_err_bg*area_err_bg);

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

    IntegratePeakWithBackground(hist, integral_min, integral_max, xmin, xmax, canvas, outFile);

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

// Calculates detector efficiency with error propagation
void efficiency_calc(double halflife, double dHalflife, double activity, double dActivity, 
                     double intensity, double dIntensity, 
                     double area, double dArea, std::ofstream &outFile, 
                     bool initial_activity, double measure_time, double decay_time)
{
    // Lambda (decay constant) and activity at measurement time
    double lambda = log(2) / halflife;
    double activity_exp;
    if(initial_activity)
        activity_exp = activity * exp(-lambda * decay_time);
    else
        activity_exp = activity * exp(lambda * decay_time);

    // I_gamma = intensity as fraction (e.g. 50% -> 0.5)
    double I_gamma = intensity / 100.0;

    // Efficiency: epsilon = peak area / (A(t) * I_gamma * t_measure)
    double efficiency = area / (activity_exp * I_gamma * measure_time);

    // ============================
    // ERROR PROPAGATION FOR EFFICIENCY
    // ============================
    // For multiplicative quantities, use relative errors.
    // In our case: epsilon ∝ area / (activity_exp * I_gamma * measure_time)
    // Log-differentiating (or using propagation rule for products/ratios):
    // (Δε/ε)^2 = (ΔA/A)^2 + (ΔA_exp/A_exp)^2 + (ΔI_gamma/I_gamma)^2 + (t_d * Δλ)^2
    // Note: the λ term appears because A_exp = A0 * exp(-λ t_d) -> logarithmic derivative gives + t_d Δλ

    // Relative error of area (from FitAndSaveToFile)
    double rel_err_area = dArea / area;

    // Relative error of initial activity (A0)
    double rel_err_activity = (dActivity > 0) ? (dActivity / activity) : 0.0;

    // Derivative of lambda with respect to T1/2: λ = ln2 / T1/2
    // dλ = (ln2 / T1/2^2) * dT1/2
    double dLambda = (log(2) / (halflife * halflife)) * dHalflife;

    // Relative error associated with decay term exp(-λ t_d) approximated as t_d * dλ
    // Explanation: A_exp = A0 * exp(-λ t_d) -> ln A_exp = ln A0 - λ t_d
    // => Δ(ln A_exp) = sqrt( (ΔA0/A0)^2 + (t_d Δλ)^2 )
    double rel_err_T12 = decay_time * dLambda;

    // Combine errors of activity (A0) and decay term
    double rel_err_activity_exp = sqrt(rel_err_activity*rel_err_activity + rel_err_T12*rel_err_T12);

    // Relative error of gamma intensity (I_gamma)
    double rel_err_Igamma   = (dIntensity / intensity);

    // Total relative error (assume independence)
    double rel_err_total = sqrt(pow(rel_err_area, 2) + pow(rel_err_activity_exp, 2) + pow(rel_err_Igamma, 2));

    // Absolute error of efficiency
    double eff_err = efficiency * rel_err_total;

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

// Analyzes peaks from a ROOT file and calculates detector efficiency
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

    std::map<int, TH1F *> histo_channels;
    std::map<int, std::vector<double>> all_resolutions;

    if (sourceData.energies.empty() || sourceData.energies.size() != sourceData.search_windows.size()) {
        std::cerr << "Error: Source data is invalid or missing (Energies/Search Windows)." << std::endl;
        input_file->Close();
        return;
    }

    std::string outputFileName = folder + "PeakAreas_Detector=" + Detector + "_Source=" + Source + ".txt";
    std::ofstream outFile(outputFileName);
    std::string pdfFileName = folder + "Fits_" + Source + "_" + Detector + ".pdf";
    TCanvas *canvas = new TCanvas("canvas", "Fits", 800, 600);
    bool firstPage = true;

    outFile << "Detector,ID,Source,Energy,Mean,Sigma,FWHM,Resolution,Area,AreaErr,"
                    << "PeakIntervalMin,PeakIntervalMax,TotalCounts,BgCounts,NetCounts,"
                    << "TotalCountsErr,BgCountsErr,NetCountsErr,"
                    << "Halflife(s),InitialActivity(Bq),dActivity,ActivityAtMeas(Bq),dActivityAtMeas,"
                    << "Intensity(%),dIntensity,MeasureTime(s),Efficiency(%),EfficiencyErr(%)"
                    << std::endl;

    for (auto i : ids)
    {
        TString name = Form("histograma_Energy_%d", i);
        histo_channels[i] = (TH1F *)input_file->Get(name);
        if (!histo_channels[i])
        {
            std::cout << "Missing histogram for det " << i << ". Continuing." << std::endl;
            continue;
        }

        std::cout << "\n>>> Detector " << i << std::endl;
        all_resolutions[i].clear();

        for (int p = 0; p < sourceData.search_windows.size(); ++p)
        {
            auto [xmin, xmax] = sourceData.search_windows[p];
            double mean_guess = sourceData.energies[p];
            double sigma_guess = sourceData.sigma_guess[p];
            double theoretical_energy = sourceData.energies[p];

            std::vector<double> area = FitAndSaveToFile(
                histo_channels[i], Detector, i, Source,
                theoretical_energy, canvas, outFile,
                all_resolutions[i], mean_guess, sigma_guess, xmin, xmax
            );

            SaveCanvasToPDF(canvas, pdfFileName, firstPage);

            efficiency_calc(
                sourceData.halflife, sourceData.dHalflife, sourceData.activity, sourceData.dActivity,
                sourceData.intensity[p], sourceData.dIntensity[p], area[0], area[1], outFile,
                sourceData.initial_activity, sourceData.measure_time, sourceData.decay_time
            );
        }
    }

    canvas->Print((pdfFileName + "]").c_str());
    input_file->Close();

    std::cout << "✅ Results have been saved in: " << outputFileName << std::endl;
}
