#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <set>
#include <algorithm>
#include <map>
#include <utility>

#include <TCanvas.h>
#include <TGraph.h>
#include <TF1.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TFile.h>
#include <TAxis.h>
#include <TFitResultPtr.h>

// ===========================================================
// Helper function to save canvas to PDF, handling multiple pages
// ===========================================================
void SaveCanvasToPDF(TCanvas *canvas, TString pdfFileName, bool &firstPage)
{
    // If this is the first page of the PDF, open with "(" to create a new PDF
    if (firstPage)
    {
        canvas->SaveAs(pdfFileName + "(");
        firstPage = false; // mark that the first page has been created
    }
    else
    {
        // For subsequent pages, just append to the existing PDF
        canvas->SaveAs(pdfFileName);
    }
}

// ===========================================================
// Main fitting function
// Reads input file with channel vs energy data for each detector
// Fits either linear or polynomial calibration curves
// Saves graphs to ROOT and PDF files
// ===========================================================
void fit(const char *fname, TString FitType, TString outputFolder)
{
    // Open input text file
    std::ifstream in(fname);
    if (!in.is_open())
    {
        std::cerr << "Error: Cannot open file " << fname << std::endl;
        return; // exit if file cannot be opened
    }

    // Create a map to store data grouped by Detector_ID
    // Key = Detector_ID, Value = vector of pairs (channel, energy)
    std::map<int, std::vector<std::pair<double, double>>> data_by_detector;
    std::string line;
    std::string source_name; // store the source name (not used for fitting here)
    int detector_id;
    double channel, energy, resolution;

    // Read the input file line by line
    while (std::getline(in, line))
    {
        std::istringstream iss(line);
        // The file format is assumed to be:
        // Source DetectorID Channel Energy Resolution
        if (iss >> source_name >> detector_id >> channel >> energy >> resolution)
        {
            // Store the channel-energy pair in the map under the corresponding detector
            data_by_detector[detector_id].push_back({channel, energy});
        }
    }
    in.close(); // close the input file

    // Check if any data was read
    if (data_by_detector.empty())
    {
        std::cerr << "Error: No data found in file " << fname << std::endl;
        return; // exit if no data
    }

    // Decide filename based on fit type
    TString filename = (FitType == "linear") ? "linear_fits" : "polynomial_fits";

    // Create a single canvas for all plots
    TCanvas *canvas = new TCanvas("canvas", filename, 800, 600);
    TString pdfFileName = outputFolder + "calibration_" + filename + ".pdf";
    bool firstPage = true; // flag to handle PDF page opening

    // Create a single ROOT file for all calibration results
    TFile *outputFile = new TFile(Form("%sEnergy_vs_Channels_All_Detectors_%s.root", outputFolder.Data(), FitType.Data()), "RECREATE");
    if (!outputFile || outputFile->IsZombie())
    {
        std::cerr << "Error: Cannot create output ROOT file." << std::endl;
        delete canvas;
        return;
    }
    outputFile->cd(); // move to output ROOT file

    // Iterate over all detectors in the map
    for (auto const &[Detector_ID, data_points] : data_by_detector)
    {
        // Extract number of data points for this detector
        int n = data_points.size();
        // Check if there are enough points to perform the fit
        if (n < 2 && FitType == "linear")
        {
            std::cerr << "Error: Not enough points for fitting for Detector " << Detector_ID << ". Skipping.\n";
            continue;
        }
        if (n < 3 && FitType == "polynomial")
        {
            std::cerr << "Error: Not enough points for polynomial fit for Detector " << Detector_ID << ". Skipping.\n";
            continue;
        }

        // Prepare vectors for x (channels) and y (energy)
        std::vector<double> x, y;
        for (const auto &point : data_points)
        {
            x.push_back(point.first);  // channel
            y.push_back(point.second); // energy
        }

        // Define the fit formula based on requested fit type
        TString fitFormula;
        if (FitType == "linear")
            fitFormula = "[0] + [1]*x"; // linear: y = a + b*x
        else if (FitType == "polynomial")
            fitFormula = "[0] + [1]*x + [2]*x^2"; // quadratic: y = a + b*x + c*x^2

        // Create TGraph and TF1 objects for plotting and fitting
        TGraph *graph = new TGraph(n, &x[0], &y[0]); // graph of data points
        TF1 *fitFunc = new TF1("fitFunc", fitFormula, x[0], x[n - 1]); // fitting function

        // Customize graph appearance
        graph->SetLineColor(kRed + 1);
        graph->SetLineWidth(3);
        graph->SetMarkerColor(kBlue + 2);
        graph->SetMarkerSize(1.5);
        graph->SetMarkerStyle(21);
        graph->SetTitle(Form("Energy vs Channels_Detector_%d", Detector_ID));
        graph->GetXaxis()->SetTitle("Channels");
        graph->GetYaxis()->SetTitle("Energy (keV)");

        // Perform the fit
        TFitResultPtr fit = graph->Fit(fitFunc, "S"); // "S" to store fit result

        // Extract fit parameters
        double intercept = fitFunc->GetParameter(0);                            // a
        double slope = fitFunc->GetParameter(1);                                // b
        double quad = (FitType == "polynomial") ? fitFunc->GetParameter(2) : 0; // c if polynomial
        double chi2ndf = (fit->Ndf() > 0) ? fit->Chi2() / fit->Ndf() : -1;      // reduced chi2

        // Draw the graph and fit on the canvas
        canvas->cd();
        graph->Draw("AP");      // A = axis, P = points
        fitFunc->Draw("SAME");  // draw fit line on same canvas

        // Add a legend with fit parameters
        TLegend *legend = new TLegend(0.15, 0.75, 0.55, 0.88);
        legend->AddEntry(graph, Form("Detector %d", Detector_ID), "p");
        legend->AddEntry((TObject *)0, Form("Slope = %.5f", slope), "");
        legend->AddEntry((TObject *)0, Form("Intercept = %.5f", intercept), "");
        legend->AddEntry((TObject *)0, Form("Quadratic = %.5f", quad), "");
        legend->AddEntry((TObject *)0, Form("#chi^{2}/NDF = %.4f", chi2ndf), "");
        legend->SetTextSize(0.03);
        legend->Draw();

        // Save the canvas to PDF (handles multiple pages)
        SaveCanvasToPDF(canvas, pdfFileName, firstPage);

        // =======================================================
        // Check if calibration entry already exists
        // =======================================================
        bool already_exists = false;
        std::ifstream in_check(Form("%scalibration_parameters.txt", outputFolder.Data()));
        if (in_check.is_open())
        {
            int existing_id;
            double s, i, d;
            std::string fit_type;
            while (in_check >> existing_id >> s >> i >> d >> fit_type)
            {
                if (existing_id == Detector_ID && fit_type == FitType)
                {
                    already_exists = true;
                    break;
                }
            }
            in_check.close();
        }

        // If entry does not exist, append it to file
        if (!already_exists)
        {
            std::ofstream out(Form("%scalibration_parameters.txt", outputFolder.Data()), std::ios::app);
            out << Detector_ID << " " << slope << " " << intercept << " " << quad << " " << FitType << std::endl;
            out.close();
        }
        else
        {
            std::cout << "⚠️ Detector " << Detector_ID << " already exists with fit type " << FitType << ". Skipping.\n";
        }

        // Save the graph, fit function, and canvas to the ROOT file
        graph->Write(Form("graph_Detector_%d", Detector_ID));
        fitFunc->Write(Form("fit_Detector_%d", Detector_ID));
        canvas->Write(Form("canvas_Detector_%d", Detector_ID));

        // Delete dynamically allocated objects to avoid memory leaks
        delete graph;
        delete fitFunc;
    }

    // Close PDF (finalize multi-page PDF)
    canvas->Clear();
    canvas->SaveAs(pdfFileName + ")");
    delete canvas;

    // Close the ROOT output file
    outputFile->Close();
    delete outputFile;
}
